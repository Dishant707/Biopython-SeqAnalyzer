<p align="center">
  <img src="https://img.shields.io/badge/Python-3.9+-blue?style=for-the-badge&logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/Biopython-1.81+-green?style=for-the-badge&logo=python&logoColor=white" alt="Biopython">
  <img src="https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge" alt="License">
  <img src="https://img.shields.io/badge/Tests-13%20Passing-success?style=for-the-badge" alt="Tests">
</p>

# 🧬 Biopython-SeqAnalyzer

A powerful Python toolkit for biological sequence analysis built on top of Biopython. Analyze FASTA files with ease using our intuitive API or command-line interface.

---

## ✨ Features

| Feature                   | Description                                                           |
| ------------------------- | --------------------------------------------------------------------- |
| 🔬 **GC Content Analysis** | Calculate GC percentage with proper handling of ambiguous nucleotides |
| 🧪 **DNA Translation**     | Translate nucleotide sequences using any NCBI codon table             |
| 🔄 **Reverse Complement**  | Generate reverse complements with IUPAC ambiguity support             |
| 🌐 **NCBI BLAST Search**   | Remote BLAST queries with automatic retry and error handling          |
| 💻 **CLI Interface**       | Full-featured command-line tool for terminal workflows                |

---

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/Dishant707/Biopython-SeqAnalyzer.git
cd Biopython-SeqAnalyzer

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```python
from src import BioSeqAnalyzer

# Load sequences
analyzer = BioSeqAnalyzer("sequences.fasta")
analyzer.load_sequences()

# Analyze GC content
gc_results = analyzer.calculate_gc_content()
for seq_id, gc_pct in gc_results.items():
    print(f"{seq_id}: {gc_pct:.2f}%")

# Translate to protein
proteins = analyzer.get_translations(table=1)

# Get reverse complements
rc_seqs = analyzer.reverse_complement_all()
```

---

## 💻 Command-Line Interface

```bash
# Show help
python run_analyzer.py --help

# Calculate GC content
python run_analyzer.py sequences.fasta --gc

# Translate sequences (bacterial codon table)
python run_analyzer.py sequences.fasta --translate --table 11

# Generate reverse complements
python run_analyzer.py sequences.fasta --reverse-complement --output rc_output.fasta

# Run BLAST search
python run_analyzer.py sequences.fasta --blast SEQ001 --database nt --hits 10

# Display sequence info
python run_analyzer.py sequences.fasta --info
```

### CLI Options

| Option                         | Description                     |
| ------------------------------ | ------------------------------- |
| `--gc`                         | Calculate GC content percentage |
| `--translate`                  | Translate to protein sequences  |
| `--reverse-complement`, `--rc` | Generate reverse complements    |
| `--blast SEQ_ID`               | Run NCBI BLAST for a sequence   |
| `--info`                       | Display sequence information    |
| `--table`, `-t`                | Codon table (default: 1)        |
| `--database`, `-d`             | BLAST database (default: nr)    |
| `--output`, `-o`               | Output file path                |

---

## 📁 Project Structure

```
Biopython-SeqAnalyzer/
├── src/
│   ├── __init__.py          # Package exports
│   └── analyzer.py          # BioSeqAnalyzer class
├── tests/
│   └── test_analyzer.py     # Unit tests (13 tests)
├── data/
│   └── test.fasta           # Sample test data
├── docs/                    # Documentation
├── run_analyzer.py          # CLI entry point
├── requirements.txt         # Dependencies
└── README.md
```

---

## 🧪 Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run with coverage
python -m pytest tests/ --cov=src --cov-report=html
```

---

## 📚 API Reference

### `BioSeqAnalyzer`

#### Constructor
```python
analyzer = BioSeqAnalyzer(file_path: str)
```

#### Methods

| Method                                    | Returns            | Description                         |
| ----------------------------------------- | ------------------ | ----------------------------------- |
| `load_sequences()`                        | `List[SeqRecord]`  | Parse FASTA file and load sequences |
| `calculate_gc_content()`                  | `Dict[str, float]` | GC percentage for each sequence     |
| `get_translations(table=1)`               | `Dict[str, str]`   | Translated protein sequences        |
| `reverse_complement_all()`                | `List[SeqRecord]`  | Reverse complement sequences        |
| `run_blast_search(seq_id, database='nr')` | `List[Dict]`       | NCBI BLAST results                  |

---

## 🔧 Requirements

- Python 3.9+
- Biopython ≥ 1.81
- pytest ≥ 7.4.0
- click ≥ 8.1.0

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

<p align="center">
  Made with ❤️ by <a href="https://github.com/Dishant707">Dishant</a>
</p>
