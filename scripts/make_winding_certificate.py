#!/usr/bin/env python3
"""Produtor do certificado intervalar de winding consumido pelo Lean.

Este programa NÃO decide nada sozinho: ele só traduz enclosures racionais
exatos para o formato que `formal/ChoptuikFormal/WindingInterval.lean`
recomputa no kernel. Toda a aritmética é feita com `fractions.Fraction`; não
existe ponto flutuante em nenhum ponto da entrada ou da saída, e qualquer
token com `.` ou `e`/`E` na posição de um número é recusado.

FORMATO DE ENTRADA (texto, uma diretiva por linha; `#` inicia comentário)
------------------------------------------------------------------------

    name <identificadorLean>   nome do `def` gerado (padrão: windingContour)
    target <inteiro>           winding esperado; se ausente, usa a soma obtida
    node-radius <racional>     raio de inflação das caixas dos nós (padrão 0)
    edge-radius <racional>     inflação extra da caixa da aresta (padrão 0)

    node <re> <im> [raio]      valor de f no nó, na ordem do contorno
    edge <reLo> <reHi> <imLo> <imHi>
                               enclosure explícito de f sobre a aresta que
                               COMEÇA no nó anterior; sobrepõe o padrão

Racionais são inteiros (`-3`) ou frações exatas (`-7/10`). O contorno fecha
implicitamente: a última aresta vai do último nó de volta ao primeiro. São
necessários pelo menos três nós.

Caixa do nó k:  [re-r, re+r] x [im-r, im+r], com r = raio do nó (ou
`node-radius`). Caixa da aresta k, quando não vem um `edge` explícito: o menor
retângulo que contém as caixas dos dois nós, inflado por `edge-radius`.

ATENÇÃO — o que este programa não pode garantir. A afirmação analítica
"`edgeBox` contém f(z) para TODO z no segmento [z_k, z_{k+1}]" é hipótese do
estimador que produziu os números; o hull inflado é apenas uma conveniência
para contornos em que o estimador já provou essa inclusão. Usar `edge-radius`
sem esse argumento é fazer ficção, não prova. Veja o cabeçalho do módulo Lean.

SAÍDAS
------

    --lean CAMINHO          bloco Lean com o `def` do certificado intervalar
                            e as regressões `by decide` (padrão: stdout)
    --increments CAMINHO    arquivo de incrementos no formato aceito por
                            scripts/verify_winding_certificate.py

Racionais não inteiros são emitidos como `mkRat n d`, e não como `n / d`:
`Rat.div` não reduz no kernel do Lean e travaria o `decide`.

O programa falha ruidosamente (código 1, mensagem no stderr) se alguma aresta
não classificar, se alguma caixa de aresta não contiver as caixas dos nós, se
alguma caixa for malformada, ou se a soma não bater com `target`.

Exemplo:

    .venv/bin/python scripts/make_winding_certificate.py \
        examples/winding-interval-square-plus1-nodes.txt \
        --increments examples/winding-interval-square-plus1-increments.txt
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_']*$")
RATIONAL = re.compile(r"^[+-]?\d+(/\d+)?$")


class CertificateError(Exception):
    """Erro fatal: o certificado não pode ser emitido."""


def parse_rational(token: str, where: str) -> Fraction:
    """Racional exato. Recusa qualquer notação de ponto flutuante."""
    if not RATIONAL.match(token):
        raise CertificateError(
            f"{where}: '{token}' não é um racional exato "
            "(use inteiro ou n/d; ponto flutuante é proibido)"
        )
    value = Fraction(token)
    if value.denominator == 0:  # inalcançável pela regex, mantido por clareza
        raise CertificateError(f"{where}: denominador nulo em '{token}'")
    return value


@dataclass(frozen=True)
class Box:
    re_lo: Fraction
    re_hi: Fraction
    im_lo: Fraction
    im_hi: Fraction

    def well_formed(self) -> bool:
        return self.re_lo <= self.re_hi and self.im_lo <= self.im_hi

    def contains(self, inner: "Box") -> bool:
        return (self.re_lo <= inner.re_lo and inner.re_hi <= self.re_hi
                and self.im_lo <= inner.im_lo and inner.im_hi <= self.im_hi)

    def strictly_pos_re(self) -> bool:
        return self.re_lo > 0

    def strictly_neg_re(self) -> bool:
        return self.re_hi < 0

    def strictly_pos_im(self) -> bool:
        return self.im_lo > 0

    def strictly_neg_im(self) -> bool:
        return self.im_hi < 0

    def inflated(self, radius: Fraction) -> "Box":
        return Box(self.re_lo - radius, self.re_hi + radius,
                   self.im_lo - radius, self.im_hi + radius)

    def __str__(self) -> str:
        return (f"[{self.re_lo}, {self.re_hi}] x "
                f"[{self.im_lo}, {self.im_hi}]")


def node_box(real: Fraction, imag: Fraction, radius: Fraction) -> Box:
    return Box(real - radius, real + radius, imag - radius, imag + radius)


def hull(left: Box, right: Box) -> Box:
    return Box(min(left.re_lo, right.re_lo), max(left.re_hi, right.re_hi),
               min(left.im_lo, right.im_lo), max(left.im_hi, right.im_hi))


def edge_crossing(start: Box, end: Box, edge: Box) -> int | None:
    """Réplica exata de `WindingEdge.edgeCrossing` do módulo Lean.

    Orientação: para o círculo percorrido no sentido anti-horário (winding
    `+1`), o cruzamento do eixo real negativo tem `Im` passando de positivo
    para negativo; logo `Im: + -> -` contribui `+1`.
    """
    if edge.strictly_pos_re():
        return 0
    if edge.strictly_pos_im() or edge.strictly_neg_im():
        return 0
    if edge.strictly_neg_re():
        # Os dois casos de MESMO sinal vêm primeiro, na mesma ordem do Lean
        # (`WindingEdge.edgeCrossing`): assim uma caixa malformada produz 0 em
        # vez de ±1, que é o erro conservador. Ver
        # `edgeCrossing_eq_zero_of_sameSignIm` em WindingInterval.lean.
        if start.strictly_pos_im() and end.strictly_pos_im():
            return 0
        if start.strictly_neg_im() and end.strictly_neg_im():
            return 0
        if start.strictly_pos_im() and end.strictly_neg_im():
            return 1
        if start.strictly_neg_im() and end.strictly_pos_im():
            return -1
        return None
    return None


@dataclass
class Contour:
    name: str
    target: int | None
    edge_radius: Fraction
    nodes: list[tuple[Fraction, Fraction, Fraction]]
    points: dict[int, tuple[Fraction, Fraction]]
    explicit_edges: dict[int, Box]


def parse_input(path: Path) -> Contour:
    name = "windingContour"
    target: int | None = None
    node_radius = Fraction(0)
    edge_radius = Fraction(0)
    nodes: list[tuple[Fraction, Fraction, Fraction]] = []
    points: dict[int, tuple[Fraction, Fraction]] = {}
    explicit_edges: dict[int, Box] = {}

    for number, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        where = f"{path}:{number}"
        parts = line.split()
        keyword, arguments = parts[0], parts[1:]

        if keyword == "name":
            if len(arguments) != 1 or not IDENTIFIER.match(arguments[0]):
                raise CertificateError(
                    f"{where}: 'name' exige um identificador Lean válido")
            name = arguments[0]
        elif keyword == "target":
            if len(arguments) != 1 or not re.match(r"^[+-]?\d+$", arguments[0]):
                raise CertificateError(f"{where}: 'target' exige um inteiro")
            target = int(arguments[0])
        elif keyword == "node-radius":
            if len(arguments) != 1:
                raise CertificateError(f"{where}: 'node-radius' exige um valor")
            node_radius = parse_rational(arguments[0], where)
            if node_radius < 0:
                raise CertificateError(f"{where}: raio negativo")
        elif keyword == "edge-radius":
            if len(arguments) != 1:
                raise CertificateError(f"{where}: 'edge-radius' exige um valor")
            edge_radius = parse_rational(arguments[0], where)
            if edge_radius < 0:
                raise CertificateError(f"{where}: raio negativo")
        elif keyword == "node":
            if len(arguments) not in (2, 3):
                raise CertificateError(
                    f"{where}: 'node' exige <re> <im> [raio]")
            real = parse_rational(arguments[0], where)
            imag = parse_rational(arguments[1], where)
            radius = (parse_rational(arguments[2], where)
                      if len(arguments) == 3 else node_radius)
            if radius < 0:
                raise CertificateError(f"{where}: raio negativo")
            nodes.append((real, imag, radius))
        elif keyword == "point":
            # Coordenada do no no plano `z`. Liga-se ao ultimo `node`
            # declarado, como `edge`. O Lean exige esses pontos desde que
            # `WindingEdge` passou a carregar `startPoint`/`endPoint`: sem eles
            # o certificado nao afirma que existe curva alguma.
            if len(arguments) != 2:
                raise CertificateError(
                    f"{where}: 'point' exige <re> <im>")
            if not nodes:
                raise CertificateError(
                    f"{where}: 'point' antes de qualquer 'node'")
            index = len(nodes) - 1
            if index in points:
                raise CertificateError(
                    f"{where}: dois pontos para o no {index}")
            points[index] = (parse_rational(arguments[0], where),
                             parse_rational(arguments[1], where))
        elif keyword == "edge":
            if len(arguments) != 4:
                raise CertificateError(
                    f"{where}: 'edge' exige <reLo> <reHi> <imLo> <imHi>")
            if not nodes:
                raise CertificateError(
                    f"{where}: 'edge' antes de qualquer 'node'")
            box = Box(*(parse_rational(token, where) for token in arguments))
            if not box.well_formed():
                raise CertificateError(f"{where}: caixa malformada {box}")
            index = len(nodes) - 1
            if index in explicit_edges:
                raise CertificateError(
                    f"{where}: duas caixas explícitas para a aresta {index}")
            explicit_edges[index] = box
        else:
            raise CertificateError(f"{where}: diretiva desconhecida '{keyword}'")

    if len(nodes) < 3:
        raise CertificateError(
            f"{path}: um contorno fechado precisa de pelo menos 3 nós "
            f"(recebidos {len(nodes)})")
    missing = [index for index in range(len(nodes)) if index not in points]
    if missing:
        raise CertificateError(
            f"{path}: sem 'point' para o(s) no(s) {missing}. Cada 'node' (que e "
            "um valor de f) precisa da sua coordenada no plano z, porque o "
            "Lean verifica a geometria do contorno, nao so o encadeamento das "
            "caixas. Acrescente 'point <re> <im>' logo apos cada 'node'.")
    return Contour(name=name, target=target, edge_radius=edge_radius,
                   nodes=nodes, points=points, explicit_edges=explicit_edges)


@dataclass
class Edge:
    start: Box
    end: Box
    box: Box
    increment: int
    start_point: tuple[Fraction, Fraction]
    end_point: tuple[Fraction, Fraction]


def build_edges(contour: Contour, edge_radius: Fraction) -> list[Edge]:
    boxes = [node_box(real, imag, radius)
             for real, imag, radius in contour.nodes]
    for index, box in enumerate(boxes):
        if not box.well_formed():
            raise CertificateError(f"nó {index}: caixa malformada {box}")

    count = len(boxes)
    edges: list[Edge] = []
    for index in range(count):
        start = boxes[index]
        end = boxes[(index + 1) % count]
        if index in contour.explicit_edges:
            box = contour.explicit_edges[index]
            source = "explícita"
        else:
            box = hull(start, end).inflated(edge_radius)
            source = "hull inflado"
        if not box.well_formed():
            raise CertificateError(f"aresta {index}: caixa malformada {box}")
        if not box.contains(start) or not box.contains(end):
            raise CertificateError(
                f"aresta {index}: a caixa da aresta ({source}) {box} não contém "
                f"as caixas dos nós {start} e {end}; o Lean rejeitaria "
                "(`shapeOk`)")
        start_point = contour.points[index]
        end_point = contour.points[(index + 1) % count]
        if start_point == end_point:
            raise CertificateError(
                f"aresta {index}: degenerada, os dois nos estao em "
                f"{start_point} no plano z; o Lean rejeitaria "
                "(`nondegenerate`)")
        increment = edge_crossing(start, end, box)
        if increment is None:
            raise CertificateError(
                f"aresta {index}: caixa {box} não classifica. Ela não está "
                "contida em {Re>0}, nem em {Im>0}, nem em {Im<0}, nem está em "
                "{Re<0} com sinais de Im estritos e opostos nos nós "
                f"({start} -> {end}). Refine o contorno ou aperte o enclosure.")
        edges.append(Edge(start=start, end=end, box=box, increment=increment,
                          start_point=start_point, end_point=end_point))
    return edges


def lean_rational(value: Fraction) -> str:
    numerator = value.numerator
    denominator = value.denominator
    if denominator == 1:
        return f"({numerator})" if numerator < 0 else str(numerator)
    head = f"({numerator})" if numerator < 0 else str(numerator)
    return f"mkRat {head} {denominator}"


def coarsen_edges(edges: list[Edge]) -> list[Edge]:
    """Une segmentos colineares orientados, preservando todos os enclosures.

    O hull contem as imagens de ambos os segmentos originais. Colinearidade
    e orientacao garantem que a uniao seja exatamente o novo segmento.
    A classificacao e recomputada, nunca inferida apenas da soma.
    """
    merged: list[Edge] = []
    for edge in edges:
        if merged:
            previous = merged[-1]
            u = tuple(b - a for a, b in zip(previous.start_point, previous.end_point))
            v = tuple(b - a for a, b in zip(edge.start_point, edge.end_point))
            collinear = u[0] * v[1] == u[1] * v[0]
            forward = u[0] * v[0] + u[1] * v[1] > 0
            if (previous.end_point == edge.start_point and previous.end == edge.start
                    and collinear and forward):
                box = hull(previous.box, edge.box)
                increment = edge_crossing(previous.start, edge.end, box)
                if increment is not None and increment == previous.increment + edge.increment:
                    merged[-1] = Edge(previous.start, edge.end, box, increment,
                                      previous.start_point, edge.end_point)
                    continue
        merged.append(edge)
    return merged if len(merged) >= 3 else edges


def lean_point(point: tuple[Fraction, Fraction]) -> str:
    return ("{ re := " + lean_rational(point[0])
            + ", im := " + lean_rational(point[1]) + " }")


def lean_box(box: Box) -> str:
    return ("{ reLo := " + lean_rational(box.re_lo)
            + ", reHi := " + lean_rational(box.re_hi)
            + ", imLo := " + lean_rational(box.im_lo)
            + ", imHi := " + lean_rational(box.im_hi) + " }")


def emit_lean(name: str, edges: list[Edge], total: int, source: Path) -> str:
    lines: list[str] = []
    lines.append("/-! Gerado por scripts/make_winding_certificate.py a partir")
    lines.append(f"de `{source}`. Não edite à mão: regenere.")
    lines.append("")
    lines.append("As três inclusões analíticas (`startBox ∋ f(z_k)`,")
    lines.append("`endBox ∋ f(z_{k+1})`, `edgeBox ∋ f` no segmento inteiro)")
    lines.append("continuam sendo hipótese do estimador que produziu os")
    lines.append("números; o que o kernel confere é tudo o mais. -/")
    lines.append(f"def {name} : Choptuik.WindingInterval.IntervalWindingCertificate where")
    lines.append("  edges :=")
    for index, edge in enumerate(edges):
        opener = "    [ " if index == 0 else "    , "
        lines.append(opener + "-- aresta " + str(index)
                     + f", incremento {edge.increment:+d}")
        lines.append("      { startPoint := " + lean_point(edge.start_point) + ",")
        lines.append("        endPoint := " + lean_point(edge.end_point) + ",")
        lines.append("        startBox := " + lean_box(edge.start) + ",")
        lines.append("        endBox := " + lean_box(edge.end) + ",")
        lines.append("        edgeBox := " + lean_box(edge.box) + " }")
    lines.append("    ]")
    lines.append(f"  target := {total}")
    lines.append("")
    lines.append(f"example : {name}.valid = true := by decide")
    lines.append(f"example : {name}.total = {total} := by decide")
    lines.append(f"example : {name}.toWinding.valid = true := by decide")
    return "\n".join(lines) + "\n"


def emit_increments(edges: list[Edge], total: int, source: Path) -> str:
    lines = [
        "# Gerado por scripts/make_winding_certificate.py a partir de",
        f"# {source}. Cada linha é o incremento de uma aresta, calculado a",
        "# partir do enclosure racional, não fornecido à mão.",
        f"# winding = {total}",
    ]
    lines.extend(str(edge.increment) for edge in edges)
    return "\n".join(lines) + "\n"


def run(arguments: argparse.Namespace) -> int:
    contour = parse_input(arguments.nodes)
    edges = build_edges(contour, contour.edge_radius)
    if arguments.coarsen:
        original_count = len(edges)
        edges = coarsen_edges(edges)
        print(f"coarsening exato: {original_count} -> {len(edges)} arestas", file=sys.stderr)
    total = sum(edge.increment for edge in edges)

    if contour.target is not None and contour.target != total:
        raise CertificateError(
            f"soma dos incrementos {total} difere do alvo declarado "
            f"{contour.target}")

    lean_text = emit_lean(contour.name, edges, total, arguments.nodes)
    if arguments.standalone:
        lean_text = ("import ChoptuikFormal.WindingInterval\n"
                     "set_option maxRecDepth 1000000\n"
                     "set_option maxHeartbeats 0\n\n" + lean_text)
    if arguments.lean is None:
        sys.stdout.write(lean_text)
    else:
        arguments.lean.write_text(lean_text)

    if arguments.increments is not None:
        arguments.increments.write_text(
            emit_increments(edges, total, arguments.nodes))

    print(f"OK: {len(edges)} arestas classificadas, winding={total}",
          file=sys.stderr)
    if not arguments.quiet:
        for index, edge in enumerate(edges):
            print(f"  aresta {index}: {edge.box} -> {edge.increment:+d}",
                  file=sys.stderr)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Emite o certificado intervalar de winding (Lean + "
                    "incrementos) a partir de enclosures racionais exatos.")
    parser.add_argument("nodes", type=Path,
                        help="arquivo de nós/arestas (ver o cabeçalho)")
    parser.add_argument("--lean", type=Path, default=None,
                        help="arquivo de saída do bloco Lean (padrão: stdout)")
    parser.add_argument("--increments", type=Path, default=None,
                        help="arquivo de incrementos para "
                             "verify_winding_certificate.py")
    parser.add_argument("--coarsen", action="store_true",
                        help="une segmentos colineares quando o hull ainda classifica")
    parser.add_argument("--standalone", action="store_true",
                        help="inclui import e limites de elaboracao para contornos grandes")
    parser.add_argument("--quiet", action="store_true", help="omite o relatorio por aresta")
    arguments = parser.parse_args()
    try:
        return run(arguments)
    except CertificateError as failure:
        print(f"ERRO: {failure}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
