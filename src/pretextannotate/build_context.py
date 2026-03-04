import argparse
from pathlib import Path


def parse_fasta_lengths(fasta_path: Path) -> list[tuple[str, int]]:
    """Return (sequence_name, sequence_length) in FASTA order."""
    records: list[tuple[str, int]] = []
    current_name: str | None = None
    current_len = 0

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records.append((current_name, current_len))
                current_name = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)

    if current_name is not None:
        records.append((current_name, current_len))

    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    return records


def parse_mapping(mapping_path: Path) -> dict[str, str]:
    """Read a 2-column TSV mapping: sequence_name<TAB>molecule_label."""
    mapping: dict[str, str] = {}
    with open(mapping_path, "r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid mapping at {mapping_path}:{line_no}. Expected 2 tab-separated columns, found {len(parts)}"
                )
            sequence_name, molecule_label = parts
            mapping[sequence_name] = molecule_label

    return mapping


def build_context_rows(
    fasta_records: list[tuple[str, int]],
    mapping: dict[str, str] | None = None,
) -> list[tuple[str, int, str]]:
    """Build rows for custom context: INSDC, length_bp, molecule_label."""
    rows: list[tuple[str, int, str]] = []
    for i, (name, length) in enumerate(fasta_records, start=1):
        molecule = mapping[name] if mapping and name in mapping else str(i)
        rows.append((name, length, molecule))
    return rows


def write_context(rows: list[tuple[str, int, str]], output_path: Path) -> None:
    with open(output_path, "w", encoding="utf-8") as out:
        out.write("# INSDC\tlength_bp\tmolecule\n")
        for name, length, molecule in rows:
            out.write(f"{name}\t{length}\t{molecule}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a custom pretextannotate context file from a FASTA file."
    )
    parser.add_argument("--fasta", required=True, type=Path, help="Input FASTA file")
    parser.add_argument(
        "--mapping",
        type=Path,
        help="Optional 2-column TSV (sequence_name<TAB>molecule_label)",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output context TSV for --sizes",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fasta_records = parse_fasta_lengths(args.fasta)
    mapping = parse_mapping(args.mapping) if args.mapping else None
    rows = build_context_rows(fasta_records, mapping)
    write_context(rows, args.output)


if __name__ == "__main__":
    main()
