"""Inspect ZIP member pickle using inert stand-ins for EVERY external class.

No external pickle-named module is imported; external functions/classes become
local data-only instances. Unlike ordinary pickle.load(), this never invokes
NetworkX, NumPy, or any callable requested by the pickle payload.
"""
from pathlib import Path
import base64
import collections
import csv
import io
import json
import pickle
import struct
import zlib

ROOT = Path(__file__).resolve().parent


class Inert:
    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        obj.args = args
        obj.kwargs = kwargs
        return obj

    def __init__(self, *args, **kwargs):
        pass

    def __setstate__(self, state):
        self.state = state


class PassiveReader(pickle.Unpickler):
    def find_class(self, module, name):
        return type("Inert_" + name, (Inert,), {"external_class": (module, name)})

    def persistent_load(self, pid):
        raise ValueError("Unsupported persistent object")


def passive_json(value):
    if isinstance(value, Inert):
        return {
            "inert_external_class": value.external_class,
            "constructor_args": passive_json(value.args),
            "constructor_kwargs": passive_json(value.kwargs),
            "state": passive_json(getattr(value, "state", None)),
        }
    if isinstance(value, type) and issubclass(value, Inert):
        return {"inert_external_class_reference": value.external_class}
    if isinstance(value, bytes):
        return {"bytes_base64": base64.b64encode(value).decode("ascii")}
    if isinstance(value, (tuple, list)):
        return [passive_json(v) for v in value]
    if isinstance(value, dict):
        return {k: passive_json(v) for k, v in value.items()}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(type(value))


member = (ROOT / "graph_annotations_member.bin").read_bytes()
assert member[:4] == b"PK\x03\x04"
header = struct.unpack("<4s5H3I2H", member[:30])
payload = zlib.decompress(member[30 + header[-2] + header[-1]:], -15)
graph = PassiveReader(io.BytesIO(payload)).load()
assert graph.external_class == ("networkx.classes.digraph", "DiGraph")
nodes = [{"node": n, **passive_json(a)} for n, a in graph.state["_node"].items()]
edges = [
    {"source": s, "target": t, **passive_json(a)}
    for s, targets in graph.state["_succ"].items() for t, a in targets.items()
]
(ROOT / "graph_from_annotations_nodes.json").write_text(json.dumps(nodes, indent=2))
(ROOT / "graph_from_annotations_edges.json").write_text(json.dumps(edges, indent=2))
for name, rows in [("nodes", nodes), ("edges", edges)]:
    with (ROOT / f"graph_from_annotations_{name}.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows({
            k: json.dumps(v) if isinstance(v, (list, dict)) else v
            for k, v in row.items()
        } for row in rows)
with (ROOT / "reactions_linex_now.csv").open() as f:
    reaction_pairs = [(r["reagent"], r["product"]) for r in csv.DictReader(f)]
graph_pairs = [(r["source"], r["target"]) for r in edges]
pair_comparison = {
    "n_graph_rows": len(graph_pairs),
    "n_reaction_rows": len(reaction_pairs),
    "ordered_pairs_identical": graph_pairs == reaction_pairs,
    "pair_multisets_identical": collections.Counter(graph_pairs) == collections.Counter(reaction_pairs),
}
(ROOT / "graph_vs_reactions_linex_pair_comparison.json").write_text(json.dumps(pair_comparison, indent=2))
print(json.dumps(pair_comparison))
