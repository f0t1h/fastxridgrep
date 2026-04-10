# fastxridgrep

A small C utility that filters FASTA/FASTQ records by read ID.

Each line in the pattern file is treated as one read ID entry; only the first
whitespace-delimited token is used.

## Build

Requirements:

- C compiler with C11 support
- `zlib`

Build:

```sh
make
```

This creates `./fastxridgrep`.

## Usage

```sh
./fastxridgrep -p ids.txt -i reads.fastq
./fastxridgrep --patterns ids.txt --input reads.fasta.gz
cat reads.fastq | ./fastxridgrep -p ids.txt
```

Options:

- `-p, --patterns` path to the read ID pattern file (required)
- `-i, --input` FASTA/FASTQ input path (`-` or omitted means stdin)

## Dependencies

Vendored under `ext/`:

- `ext/kseq.h` from klib
- `ext/aho-corasick` from mischasan/aho-corasick
