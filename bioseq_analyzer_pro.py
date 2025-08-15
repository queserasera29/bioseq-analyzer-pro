# BioSeq Analyzer Pro – PyQt5
# DNA→Protein, ORFs, Live Mutation Sandbox, Hydropathy Plot, CRISPR PAM Finder, Primer Tm, FASTA/Report Export
# Author: Assistant for Ranjana
# Python: 3.9+ (tested up to 3.12)
# Dependencies: pyqt5, matplotlib (installs numpy)

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QTextEdit, QLineEdit, QTabWidget, QFileDialog,
    QSpinBox, QComboBox, QTableWidget, QTableWidgetItem, QMessageBox, QSplitter,
    QGroupBox, QAction
, QDialog, QTextBrowser)

import sys, os
from collections import Counter

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

APP_TITLE = "BioSeq Analyzer Pro – PyQt5"

# -------------------------
# Core Bioinformatics Utils
# -------------------------
DNA_ALPHABET = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

HYDROPATHY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4,
    'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
    'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5, '*': 0.0, 'X': 0.0
}

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def clean_dna(seq: str) -> str:
    s = ''.join([c for c in seq if c.upper() in 'ACGT'])
    return s.upper()

def reverse_complement(dna: str) -> str:
    return dna.translate(COMPLEMENT)[::-1].upper()

def transcribe(dna: str) -> str:
    return dna.replace('T', 'U')

def translate(dna: str, frame: int = 0) -> str:
    dna = dna[frame:]
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)

def gc_content(dna: str) -> float:
    if not dna:
        return 0.0
    gc = dna.count('G') + dna.count('C')
    return (gc / len(dna)) * 100.0

def find_motif_positions(seq: str, motif: str) -> list:
    pos = []
    if not motif:
        return pos
    i = 0
    m = motif.upper()
    s = seq.upper()
    while True:
        j = s.find(m, i)
        if j == -1:
            break
        pos.append(j)
        i = j + 1
    return pos

def codon_usage(dna: str) -> dict:
    counts = Counter()
    for i in range(0, len(dna) - 2, 3):
        counts[dna[i:i+3]] += 1
    total = sum(counts.values()) or 1
    return {c: counts[c] / total for c in sorted(counts)}

def find_orfs(dna: str, min_len_aa: int = 30) -> list:
    results = []
    for frame in range(3):
        i = frame
        while i <= len(dna) - 3:
            if dna[i:i+3] == START_CODON:
                j = i
                while j <= len(dna) - 3:
                    c = dna[j:j+3]
                    if c in STOP_CODONS:
                        orf_dna = dna[i:j+3]
                        aa = translate(orf_dna)
                        aa = aa.split('*')[0]
                        if len(aa) >= min_len_aa:
                            results.append({
                                'frame': frame,
                                'start_nt': i + 1,
                                'end_nt': j + 3,
                                'length_aa': len(aa),
                                'protein': aa,
                                'dna': orf_dna
                            })
                        i = j + 3
                        break
                    j += 3
            i += 3
    return sorted(results, key=lambda x: (-x['length_aa'], x['start_nt']))

def hydropathy_profile(protein: str, window: int = 9) -> list:
    if window < 1:
        window = 1
    vals = [HYDROPATHY.get(a, 0.0) for a in protein]
    n = len(vals)
    if n == 0:
        return []
    half = window // 2
    scores = []
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        w = max(1, b - a)
        scores.append(sum(vals[a:b]) / w)
    return scores

def find_pam_sites(dna: str, pam: str = 'NGG') -> list:
    dna = dna.upper()
    rcdna = reverse_complement(dna)

    def match(seq, motif):
        for a, b in zip(seq, motif):
            if b == 'N':
                continue
            if a != b:
                return False
        return True

    sites = []
    # + strand
    for i in range(len(dna) - len(pam) + 1):
        if match(dna[i:i+len(pam)], pam):
            guide_start = max(0, i - 20)
            guide = dna[guide_start:i]
            sites.append({'strand': '+', 'pam_start': i+1, 'pam': dna[i:i+3], 'guide': guide})
    # - strand (rev PAM on forward)
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    rev_pam = ''.join(comp.get(b, 'N') for b in pam)[::-1]
    for i in range(len(dna) - len(rev_pam) + 1):
        if match(dna[i:i+len(rev_pam)], rev_pam):
            guide = reverse_complement(dna[i+3:i+3+20])
            sites.append({'strand': '-', 'pam_start': i+1, 'pam': dna[i:i+3], 'guide': guide})
    return sites

def wrap_text(seq: str, width: int = 60) -> str:
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

# -------------------------
# Matplotlib Canvas Widget
# -------------------------
class MplCanvas(FigureCanvas):
    def __init__(self, width=5, height=3, dpi=100, parent=None):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        if parent:
            self.setParent(parent)

# -------------------------
# Main Window
# -------------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_TITLE)
        self.resize(1200, 800)
        self.original_dna = ''
        self.current_dna = ''

        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        self._build_menu()
        self._build_tab_analyze()
        self._build_tab_mutate()
        self._build_tab_crispr_primer()

        self.statusBar().showMessage('Ready')

    # ---- Menu ----
    def _build_menu(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu('&File')

        act_open = QAction('Open Sequence…', self)
        act_open.triggered.connect(self.open_sequence)
        file_menu.addAction(act_open)

        act_save_report = QAction('Save Report (.txt)…', self)
        act_save_report.triggered.connect(self.save_report)
        file_menu.addAction(act_save_report)

        act_export_fasta = QAction('Export FASTA (original/mutated)…', self)
        act_export_fasta.triggered.connect(self.export_fasta)
        file_menu.addAction(act_export_fasta)

        file_menu.addSeparator()
        act_quit = QAction('Quit', self)
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

        help_menu = menubar.addMenu('&Help')
        act_tutorial = QAction('Quick Tutorial', self)
        act_tutorial.triggered.connect(self.show_tutorial)
        help_menu.addAction(act_tutorial)

        act_about = QAction('About', self)
        act_about.triggered.connect(self.show_about)
        help_menu.addAction(act_about)

    # ---- Tab 1: Analyze ----
    def _build_tab_analyze(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # Controls row
        controls = QHBoxLayout()
        self.txt_dna = QTextEdit()
        self.txt_dna.setPlaceholderText('Paste or type DNA sequence (A/C/G/T). Headers or whitespace allowed; non-ACGT will be removed.')
        self.txt_dna.setFixedHeight(120)

        buttons = QVBoxLayout()
        btn_load = QPushButton('Load .txt/.fa')
        btn_load.clicked.connect(self.open_sequence)
        btn_demo = QPushButton('Use Demo')
        btn_demo.clicked.connect(self.load_demo)
        btn_clear = QPushButton('Clear')
        btn_clear.clicked.connect(lambda: self.txt_dna.clear())
        btn_analyze = QPushButton('Analyze')
        btn_analyze.clicked.connect(self.analyze)
        for b in (btn_load, btn_demo, btn_clear, btn_analyze):
            b.setMinimumWidth(140)
            buttons.addWidget(b)
        buttons.addStretch(1)

        controls_left = QVBoxLayout()
        controls_left.addWidget(QLabel('Input DNA Sequence'))
        controls_left.addWidget(self.txt_dna)

        controls_layout = QHBoxLayout()
        controls_layout.addLayout(controls_left, 4)
        controls_layout.addLayout(buttons, 1)

        layout.addLayout(controls_layout)

        # Frame + Motif
        fm = QHBoxLayout()
        fm.addWidget(QLabel('Reading frame:'))
        self.spin_frame = QSpinBox()
        self.spin_frame.setRange(0,2)
        fm.addWidget(self.spin_frame)
        fm.addSpacing(20)
        fm.addWidget(QLabel('Motif:'))
        self.edit_motif = QLineEdit()
        self.edit_motif.setPlaceholderText('e.g., ATG')
        fm.addWidget(self.edit_motif)
        layout.addLayout(fm)

        # Splitter for outputs
        splitter = QSplitter(Qt.Vertical)

        # Overview
        box_over = QGroupBox('Overview')
        over_layout = QVBoxLayout(box_over)
        self.txt_overview = QTextEdit()
        self.txt_overview.setReadOnly(True)
        over_layout.addWidget(self.txt_overview)
        splitter.addWidget(box_over)

        # Protein + ORFs + Hydropathy
        lower = QWidget()
        lower_lay = QGridLayout(lower)

        self.txt_protein = QTextEdit()
        self.txt_protein.setReadOnly(True)
        lower_lay.addWidget(QLabel('Protein (stop codons highlighted)'), 0, 0)
        lower_lay.addWidget(self.txt_protein, 1, 0)

        self.tbl_orf = QTableWidget(0, 4)
        self.tbl_orf.setHorizontalHeaderLabels(['Frame','Start (nt)','End (nt)','Length (aa)'])
        self.tbl_orf.horizontalHeader().setStretchLastSection(True)
        lower_lay.addWidget(QLabel('ORFs (≥30 aa)'), 0, 1)
        lower_lay.addWidget(self.tbl_orf, 1, 1)

        self.canvas = MplCanvas(width=5, height=2.6, dpi=100)
        lower_lay.addWidget(QLabel('Hydropathy (Kyte-Doolittle, main ORF)'), 2, 0, 1, 2)
        lower_lay.addWidget(self.canvas, 3, 0, 1, 2)

        splitter.addWidget(lower)

        layout.addWidget(splitter)

        self.tabs.addTab(tab, 'Sequence Analysis')

    # ---- Tab 2: Mutation Sandbox ----
    def _build_tab_mutate(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        row = QHBoxLayout()
        row.addWidget(QLabel('Position (1-based nt):'))
        self.spin_pos = QSpinBox()
        self.spin_pos.setRange(1, 10000000)
        row.addWidget(self.spin_pos)
        row.addSpacing(10)
        row.addWidget(QLabel('New base:'))
        self.combo_base = QComboBox()
        self.combo_base.addItems(['A','C','G','T'])
        row.addWidget(self.combo_base)

        btn_apply = QPushButton('Apply Mutation')
        btn_apply.clicked.connect(self.apply_mutation)
        row.addWidget(btn_apply)

        btn_reset = QPushButton('Reset')
        btn_reset.clicked.connect(self.reset_mutation)
        row.addWidget(btn_reset)
        row.addStretch(1)
        layout.addLayout(row)

        self.txt_mut_protein = QTextEdit()
        self.txt_mut_protein.setReadOnly(True)
        layout.addWidget(QLabel('Mutated Protein (differences highlighted)'))
        layout.addWidget(self.txt_mut_protein)

        self.canvas_mut = MplCanvas(width=6, height=2.6, dpi=100)
        layout.addWidget(QLabel('Hydropathy (mutated main ORF)'))
        layout.addWidget(self.canvas_mut)

        self.tabs.addTab(tab, 'Mutation Sandbox')

    # ---- Tab 3: CRISPR + Primer ----
    def _build_tab_crispr_primer(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # CRISPR
        crispr_box = QGroupBox('CRISPR PAM Finder (SpCas9 default NGG)')
        cr_lay = QHBoxLayout(crispr_box)
        cr_lay.addWidget(QLabel('PAM motif:'))
        self.edit_pam = QLineEdit('NGG')
        self.edit_pam.setFixedWidth(80)
        cr_lay.addWidget(self.edit_pam)
        btn_find = QPushButton('Find PAM Sites')
        btn_find.clicked.connect(self.find_pams)
        cr_lay.addWidget(btn_find)
        cr_lay.addStretch(1)

        self.tbl_pam = QTableWidget(0, 4)
        self.tbl_pam.setHorizontalHeaderLabels(['Strand','PAM Start','PAM','20nt Guide'])
        self.tbl_pam.horizontalHeader().setStretchLastSection(True)

        # Primer
        primer_box = QGroupBox('Primer Tm Calculator')
        pr_lay = QHBoxLayout(primer_box)
        pr_lay.addWidget(QLabel('Primer sequence:'))
        self.edit_primer = QLineEdit()
        self.edit_primer.setPlaceholderText('e.g., ATGCGTACGTTAG...')
        pr_lay.addWidget(self.edit_primer)
        btn_tm = QPushButton('Compute Tm')
        btn_tm.clicked.connect(self.compute_tm)
        pr_lay.addWidget(btn_tm)
        self.lbl_tm = QLabel('')
        pr_lay.addWidget(self.lbl_tm)
        pr_lay.addStretch(1)

        layout.addWidget(crispr_box)
        layout.addWidget(self.tbl_pam)
        layout.addWidget(primer_box)

        self.tabs.addTab(tab, 'CRISPR & Primer')

    # ---- Actions ----
    def open_sequence(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Open sequence file', '', 'Sequence files (*.txt *.fa *.fasta);;All files (*)')
        if not path:
            return
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                raw = f.read()
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not open file:\n{e}')
            return
        if raw.startswith('>'):
            lines = [ln.strip() for ln in raw.splitlines() if not ln.startswith('>')]
            raw = ''.join(lines)
        self.txt_dna.setText(raw.strip())

    def load_demo(self):
        demo = (
            "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
            "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGAT"
            "GTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTT"
        )
        self.txt_dna.setText(demo)

    def analyze(self):
        raw = self.txt_dna.toPlainText().strip()
        dna = clean_dna(raw)
        self.original_dna = dna
        self.current_dna = dna
        if not dna:
            QMessageBox.information(self, 'Info', 'Please enter a DNA sequence (A/C/G/T).')
            return
        frame = self.spin_frame.value()
        motif = self.edit_motif.text().strip().upper()

        rcomp = reverse_complement(dna)
        rna = transcribe(dna)
        prot = translate(dna, frame)
        gc = gc_content(dna)
        motif_pos = find_motif_positions(dna, motif) if motif else []
        orfs = find_orfs(dna)
        cu = codon_usage(dna)

        # Overview
        over = []
        over.append(f"Length (nt): {len(dna)}")
        over.append(f"GC%: {gc:.2f}")
        if motif:
            over.append(f"Motif '{motif}': {len(motif_pos)} at positions {[p+1 for p in motif_pos]}")
        over.append("")
        over.append("DNA (cleaned):\n" + wrap_text(dna))
        over.append("")
        over.append("Reverse Complement:\n" + wrap_text(rcomp))
        over.append("")
        over.append(f"mRNA (frame {frame}):\n" + wrap_text(rna))
        over.append("")
        # Codon usage brief (top few)
        if cu:
            top = sorted(cu.items(), key=lambda x: -x[1])[:8]
            over.append("Top codon usage (fraction): " + ', '.join([f"{k}:{v:.2f}" for k,v in top]))
        self.txt_overview.setPlainText('\n'.join(over))

        # Protein with stop highlighted
        self.render_protein(self.txt_protein, prot)

        # ORF table
        self.tbl_orf.setRowCount(0)
        for it in orfs:
            r = self.tbl_orf.rowCount()
            self.tbl_orf.insertRow(r)
            self.tbl_orf.setItem(r, 0, QTableWidgetItem(str(it['frame'])))
            self.tbl_orf.setItem(r, 1, QTableWidgetItem(str(it['start_nt'])))
            self.tbl_orf.setItem(r, 2, QTableWidgetItem(str(it['end_nt'])))
            self.tbl_orf.setItem(r, 3, QTableWidgetItem(str(it['length_aa'])))

        # Hydropathy main ORF (up to first stop)
        main_prot = prot.split('*')[0]
        self.plot_hydropathy(self.canvas, main_prot)
        self.statusBar().showMessage('Analysis complete')

        # Reset mutation tab displays
        self.txt_mut_protein.clear()
        self.canvas_mut.ax.clear(); self.canvas_mut.ax.set_title(''); self.canvas_mut.draw()

        # Clear PAM table until requested
        self.tbl_pam.setRowCount(0)

    def render_protein(self, widget: QTextEdit, protein: str, diffs: set = None):
        widget.clear()
        cursor = widget.textCursor()
        fmt_default = QtGui.QTextCharFormat()
        fmt_stop = QtGui.QTextCharFormat(); fmt_stop.setForeground(QtGui.QBrush(Qt.red))
        fmt_diff = QtGui.QTextCharFormat(); fmt_diff.setBackground(QtGui.QBrush(QtGui.QColor('#fff4cc')))

        for i, ch in enumerate(protein):
            if diffs and i in diffs:
                cursor.setCharFormat(fmt_diff)
            elif ch == '*':
                cursor.setCharFormat(fmt_stop)
            else:
                cursor.setCharFormat(fmt_default)
            cursor.insertText(ch)
        cursor.insertText('\n')

    def plot_hydropathy(self, canvas: MplCanvas, protein: str):
        canvas.ax.clear()
        canvas.ax.set_xlabel('Position (aa)')
        canvas.ax.set_ylabel('Hydropathy (Kyte-Doolittle)')
        if protein:
            scores = hydropathy_profile(protein, window=9)
            xs = list(range(1, len(scores)+1))
            canvas.ax.plot(xs, scores)
            canvas.ax.set_title(f'Len={len(protein)} aa')
        else:
            canvas.ax.set_title('No protein')
        canvas.draw()

    def apply_mutation(self):
        if not self.original_dna:
            QMessageBox.information(self, 'Info', 'Run Analyze first.')
            return
        pos = self.spin_pos.value()
        base = self.combo_base.currentText()
        if pos < 1 or pos > len(self.original_dna):
            QMessageBox.critical(self, 'Error', 'Position out of range.')
            return
        dna_list = list(self.original_dna)
        dna_list[pos-1] = base
        self.current_dna = ''.join(dna_list)

        frame = self.spin_frame.value()
        prot_ref = translate(self.original_dna, frame)
        prot_mut = translate(self.current_dna, frame)

        diffs = set(i for i,(a,b) in enumerate(zip(prot_ref, prot_mut)) if a != b)
        self.render_protein(self.txt_mut_protein, prot_mut, diffs)
        self.plot_hydropathy(self.canvas_mut, prot_mut.split('*')[0])
        self.statusBar().showMessage(f'Mutation applied at nt {pos} → {base}')

    def reset_mutation(self):
        if not self.original_dna:
            return
        self.current_dna = self.original_dna
        frame = self.spin_frame.value()
        prot = translate(self.original_dna, frame)
        self.render_protein(self.txt_mut_protein, prot)
        self.plot_hydropathy(self.canvas_mut, prot.split('*')[0])
        self.statusBar().showMessage('Reset to original sequence')

    def find_pams(self):
        dna = clean_dna(self.txt_dna.toPlainText())
        if not dna:
            QMessageBox.information(self, 'Info', 'Enter a DNA sequence on the Analysis tab.')
            return
        pam = self.edit_pam.text().strip().upper() or 'NGG'
        sites = find_pam_sites(dna, pam=pam)
        self.tbl_pam.setRowCount(0)
        for s in sites:
            r = self.tbl_pam.rowCount(); self.tbl_pam.insertRow(r)
            self.tbl_pam.setItem(r, 0, QTableWidgetItem(s['strand']))
            self.tbl_pam.setItem(r, 1, QTableWidgetItem(str(s['pam_start'])))
            self.tbl_pam.setItem(r, 2, QTableWidgetItem(s['pam']))
            self.tbl_pam.setItem(r, 3, QTableWidgetItem(s['guide']))
        if not sites:
            QMessageBox.information(self, 'Result', 'No PAM sites found with the chosen motif.')

    def compute_tm(self):
        primer = self.edit_primer.text().strip().upper()
        if not primer:
            QMessageBox.information(self, 'Info', 'Enter a primer sequence.')
            return
        gc = gc_content(primer)
        # Wallace rule for <14, otherwise SantaLucia-like approximation
        if len(primer) < 14:
            tm = (primer.count('A') + primer.count('T'))*2 + (primer.count('G') + primer.count('C'))*4
        else:
            tm = 64.9 + 41 * (primer.count('G') + primer.count('C') - 16.4) / len(primer)
        self.lbl_tm.setText(f" Len: {len(primer)} | GC%: {gc:.2f} | Tm: {tm:.2f}°C")

    def save_report(self):
        if not (self.original_dna or self.current_dna):
            QMessageBox.information(self, 'Info', 'Analyze a sequence first.')
            return
        path, _ = QFileDialog.getSaveFileName(self, 'Save report', 'bioseq_report.txt', 'Text files (*.txt)')
        if not path:
            return
        dna = self.current_dna or self.original_dna
        frame = self.spin_frame.value()
        report = []
        report.append(APP_TITLE)
        report.append('='*len(APP_TITLE))
        report.append('')
        report.append(f'Length (nt): {len(dna)}')
        report.append(f'GC%: {gc_content(dna):.2f}')
        report.append('DNA (cleaned):')
        report.append(wrap_text(dna))
        report.append('')
        report.append(f'Protein (frame {frame}):')
        prot = translate(dna, frame)
        report.append(prot)
        report.append('')
        report.append('Top ORFs (≥30 aa):')
        for it in find_orfs(dna):
            report.append(f"  Frame {it['frame']} | {it['start_nt']}-{it['end_nt']} | {it['length_aa']} aa")
        report.append('')
        report.append('CRISPR PAM (NGG by default):')
        for s in find_pam_sites(dna):
            report.append(f"  {s['strand']} strand | PAM@{s['pam_start']} {s['pam']} | guide: {s['guide']}")
        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write('\n'.join(report))
            QMessageBox.information(self, 'Saved', f'Report saved to\n{path}')
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not save report:\n{e}')

    def export_fasta(self):
        if not (self.original_dna or self.current_dna):
            QMessageBox.information(self, 'Info', 'Analyze a sequence first.')
            return
        dna = self.current_dna or self.original_dna
        base, _ = QFileDialog.getSaveFileName(self, 'Export FASTA', 'sequence.fasta', 'FASTA (*.fasta *.fa)')
        if not base:
            return
        try:
            with open(base, 'w', encoding='utf-8') as f:
                f.write('>bioseq_original_or_current\n')
                for line in wrap_text(dna).split('\n'):
                    f.write(line+'\n')
            QMessageBox.information(self, 'Saved', f'FASTA exported to\n{base}')
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not export FASTA:\n{e}')

    def show_tutorial(self):
        html = """
        <style>
          body { font-family: -apple-system, Segoe UI, Roboto, Arial, sans-serif; line-height: 1.35; }
          code, pre { background:#f6f8fa; padding:2px 4px; border-radius:4px; }
          h2{margin-top:0}
          li{margin:4px 0}
        </style>
        <h2>Quick Tutorial</h2>
        <ol>
          <li><b>Sequence Analysis:</b> Paste sequence (A/C/G/T) or click <i>Use Demo</i> → <i>Analyze</i>. View Overview, Protein (stops in red), ORFs, and Hydropathy.</li>
          <li><b>Mutation Sandbox:</b> Enter 1-based position + new base → <i>Apply Mutation</i>. Differing amino acids are highlighted. Use <i>Reset</i> to revert.</li>
          <li><b>CRISPR & Primer:</b> Set PAM (default NGG) → <i>Find PAM Sites</i>. Enter primer → <i>Compute Tm</i> to see length, GC%, and Tm.</li>
          <li><b>File → Save Report / Export FASTA</b> for submissions.</li>
        </ol>
        <h3>Result Interpretation (Cheat‑Sheet)</h3>
        <ul>
          <li><b>GC%:</b> Higher GC → more stable DNA; affects PCR/Tm.</li>
          <li><b>ORFs:</b> Long ORF in one frame suggests coding region; short ORFs (&lt;50 aa) may be false positives.</li>
          <li><b>Hydropathy:</b> Positive peaks → hydrophobic (possible transmembrane); negative → hydrophilic/soluble.</li>
          <li><b>Mutations:</b> No AA change = synonymous; AA change = missense; "*" gained = nonsense (truncation).</li>
          <li><b>CRISPR:</b> NGG required for SpCas9; pick guides with ~40–60% GC and validate off‑targets externally.</li>
          <li><b>Primers:</b> 18–24 nt, 40–60% GC, GC clamp at 3′; avoid strong self‑complementarity.</li>
        </ul>
        <p><b>Report tip:</b> Include sequence length & GC%, main ORF coordinates, mutation effect, hydropathy figure, chosen CRISPR guides, and primer Tm.</p>
        """
        dlg = QDialog(self)
        dlg.setWindowTitle('Quick Tutorial')
        dlg.resize(760, 520)
        lay = QVBoxLayout(dlg)
        tb = QTextBrowser(); tb.setHtml(html)
        lay.addWidget(tb)
        btn = QPushButton('Close'); btn.clicked.connect(dlg.accept)
        lay.addWidget(btn, alignment=Qt.AlignRight)
        dlg.exec_()

    def show_about(self):
        QMessageBox.information(self, 'About', (
            f"{APP_TITLE}\n\n"
            "Features:\n"
            "• DNA cleanup, reverse complement, mRNA\n"
            "• Translation (frame 0–2), stop codons highlighted\n"
            "• ORF finder (≥30 aa) with table\n"
            "• Hydropathy plots (Kyte–Doolittle)\n"
            "• Live Mutation Sandbox with AA diffs\n"
            "• CRISPR PAM (NGG default) with 20nt guides\n"
            "• Primer Tm calculator\n"
            "• Save text report, export FASTA\n\n"
            "Dependencies: PyQt5, matplotlib"
        ))

# -------------------------
# Main Entrypoint
# -------------------------
if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())
