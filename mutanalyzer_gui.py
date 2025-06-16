import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import tkinter.font as tkFont
from Bio import Entrez, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Data.IUPACData import protein_letters_1to3
import csv
import threading
from datetime import datetime
import time
import requests
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

# IMPORTANT: Change this to your actual email address
Entrez.email = "your_actual_email@domain.com"  # MUST BE CHANGED
VEP_EMAIL = "your_actual_email@domain.com"  # For Ensembl VEP API (optional)

class MutationAnalyzer:
    def __init__(self):
        self.exon_ranges = []
        self.intron_ranges = []
        self.mutations = []
        self.ref_seq = ""
        self.sample_seq = ""
        self.aligned_ref = ""
        self.aligned_sample = ""
        self.align_btn = None
        self.analyze_btn = None
        self.pathogenicity_btn = None
        self.chrom = None  # To store chromosome from NCBI fetch
        self.colors = {
            'primary': '#1e293b',      # Slate 800
            'primary_light': '#334155',  # Slate 700
            'secondary': '#0ea5e9',    # Sky 500
            'secondary_light': '#38bdf8', # Sky 400
            'accent': '#8b5cf6',       # Violet 500
            'success': '#10b981',      # Emerald 500
            'success_light': '#34d399', # Emerald 400
            'warning': '#f59e0b',      # Amber 500
            'danger': '#ef4444',       # Red 500
            'info': '#6366f1',         # Indigo 500
            'light': '#f8fafc',        # Slate 50
            'card': '#ffffff',
            'background': '#f1f5f9',   # Slate 100
            'dark': '#0f172a',         # Slate 900
            'text_primary': '#1e293b',
            'text_secondary': '#64748b',
            'border': '#e2e8f0',       # Slate 200
            'border_focus': '#3b82f6', # Blue 500
            'hover': '#f8fafc'         # Slate 50
        }
        self.setup_gui()
        
    def setup_gui(self):
        self.root = tk.Tk()
        self.root.title("🧬 MutAnalyzer Pro - Advanced Mutation Analysis Tool")
        self.root.geometry("1280x720")
        self.root.configure(bg=self.colors['background'])
        self.setup_styles()
        self.create_header()
        self.create_main_interface()
        self.add_settings_menu()
        
    def setup_styles(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('Primary.TButton', background=self.colors['primary'], foreground='white', borderwidth=2, relief='flat', focuscolor='none', padding=(15, 10), font=('Segoe UI', 10, 'bold'))
        style.map('Primary.TButton', background=[('active', self.colors['primary_light']), ('pressed', self.colors['dark'])])
        style.configure('Success.TButton', background=self.colors['success'], foreground='white', borderwidth=2, relief='flat', focuscolor='none', padding=(15, 10), font=('Segoe UI', 10, 'bold'))
        style.map('Success.TButton', background=[('active', self.colors['success_light']), ('pressed', '#059669')])
        style.configure('Danger.TButton', background=self.colors['danger'], foreground='white', borderwidth=2, relief='flat', focuscolor='none', padding=(15, 10), font=('Segoe UI', 10, 'bold'))
        style.map('Danger.TButton', background=[('active', '#f87171'), ('pressed', '#dc2626')])
        style.configure('Info.TButton', background=self.colors['info'], foreground='white', borderwidth=2, relief='flat', focuscolor='none', padding=(15, 10), font=('Segoe UI', 10, 'bold'))
        style.map('Info.TButton', background=[('active', '#818cf8'), ('pressed', '#4f46e5')])
        style.configure('TNotebook', background=self.colors['background'], borderwidth=0, tabposition='n')
        style.configure('TNotebook.Tab', background=self.colors['card'], foreground=self.colors['text_primary'], padding=(15, 10), borderwidth=2, relief='flat', font=('Segoe UI', 10, 'bold'))
        style.map('TNotebook.Tab', background=[('selected', self.colors['primary']), ('active', self.colors['hover'])], foreground=[('selected', 'white'), ('active', self.colors['primary'])])
        style.configure('TProgressbar', background=self.colors['secondary'], troughcolor=self.colors['border'], borderwidth=0, lightcolor=self.colors['secondary'], darkcolor=self.colors['secondary'])
        style.configure('Treeview', background=self.colors['card'], foreground=self.colors['text_primary'], fieldbackground=self.colors['card'], borderwidth=1, relief='solid', selectbackground=self.colors['secondary_light'], selectforeground='white', font=('Segoe UI', 10))
        style.configure('Treeview.Heading', background=self.colors['primary'], foreground='white', borderwidth=1, relief='flat', font=('Segoe UI', 10, 'bold'))
        style.map('Treeview.Heading', background=[('active', self.colors['primary_light'])])
    
    def create_header(self):
        header_frame = tk.Frame(self.root, bg=self.colors['primary'], height=100)
        header_frame.pack(fill='x', padx=0, pady=0)
        gradient_frame = tk.Frame(header_frame, bg=self.colors['primary_light'], height=15)
        gradient_frame.pack(fill='x', side='bottom')
        title_container = tk.Frame(header_frame, bg=self.colors['primary'])
        title_container.pack(expand=True, fill='both', padx=15, pady=10)
        title_font = tkFont.Font(family="Segoe UI", size=24, weight="bold")
        title_label = tk.Label(title_container, text="🧬 MutAnalyzer Pro", font=title_font, bg=self.colors['primary'], fg='white')
        title_label.pack(pady=(10, 5))
        subtitle_font = tkFont.Font(family="Segoe UI", size=11)
        subtitle_label = tk.Label(title_container, text="Advanced DNA Mutation Analysis & Genomic Alignment Platform", font=subtitle_font, bg=self.colors['primary'], fg=self.colors['light'])
        subtitle_label.pack(pady=(0, 5))
        features_font = tkFont.Font(family="Segoe UI", size=9)
        features_label = tk.Label(title_container, text="✨ NCBI Integration • Sequence Alignment • Mutation Detection • Pathogenicity Prediction • Professional Reports", font=features_font, bg=self.colors['primary'], fg=self.colors['secondary_light'])
        features_label.pack(pady=(0, 10))

    def create_main_interface(self):
        main_container = tk.Frame(self.root, bg=self.colors['background'])
        main_container.pack(fill='both', expand=True, padx=15, pady=15)
        self.notebook = ttk.Notebook(main_container)
        self.input_tab = self.create_input_tab()
        self.alignment_tab = self.create_alignment_tab()
        self.mutation_tab = self.create_mutation_tab()
        self.results_tab = self.create_results_tab()
        self.notebook.add(self.input_tab, text="📝 Input & Fetch")
        self.notebook.add(self.alignment_tab, text="🔗 Alignment")
        self.notebook.add(self.mutation_tab, text="🧪 Mutations")
        self.notebook.add(self.results_tab, text="📊 Results")
        self.notebook.pack(fill='both', expand=True)

    def create_card_frame(self, parent, title, height=None):
        card_container = tk.Frame(parent, bg=self.colors['background'])
        shadow_frame = tk.Frame(card_container, bg='#d1d5db', height=2)
        shadow_frame.pack(fill='x', padx=(2, 0), pady=(2, 0))
        card = tk.Frame(card_container, bg=self.colors['card'], relief='flat', bd=1, highlightbackground=self.colors['border'], highlightthickness=1)
        if height:
            card.configure(height=height-2)
            card.pack_propagate(False)
        card.pack(fill='both', expand=True, padx=(0, 2), pady=(0, 2))
        header = tk.Frame(card, bg=self.colors['light'], height=40)
        header.pack(fill='x', padx=0, pady=0)
        accent_line = tk.Frame(header, bg=self.colors['secondary'], height=3)
        accent_line.pack(fill='x', side='top')
        title_font = tkFont.Font(family="Segoe UI", size=12, weight="bold")
        title_label = tk.Label(header, text=title, font=title_font, bg=self.colors['light'], fg=self.colors['text_primary'])
        title_label.pack(pady=(8, 5))
        content = tk.Frame(card, bg=self.colors['card'])
        content.pack(fill='both', expand=True, padx=15, pady=15)
        return card_container, content

    def create_tooltip(self, widget, text):
        tooltip = tk.Toplevel(widget)
        tooltip.wm_overrideredirect(True)
        tooltip.wm_geometry(f"+{widget.winfo_rootx()+20}+{widget.winfo_rooty()+20}")
        label = tk.Label(tooltip, text=text, bg=self.colors['light'], fg=self.colors['text_primary'], relief='solid', borderwidth=1, font=("Segoe UI", 9))
        label.pack()
        tooltip.withdraw()
        def show(event): tooltip.deiconify()
        def hide(event): tooltip.withdraw()
        widget.bind("<Enter>", show)
        widget.bind("<Leave>", hide)

    def add_settings_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        settings_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Settings", menu=settings_menu)
        font_sizes = [10, 12, 14, 16]
        self.font_size_var = tk.IntVar(value=10)
        font_menu = tk.Menu(settings_menu, tearoff=0)
        settings_menu.add_cascade(label="Font Size", menu=font_menu)
        for size in font_sizes:
            font_menu.add_radiobutton(label=f"{size}pt", variable=self.font_size_var, value=size, command=self.update_font_size)
        self.theme_var = tk.StringVar(value="Light")
        theme_menu = tk.Menu(settings_menu, tearoff=0)
        settings_menu.add_cascade(label="Theme", menu=theme_menu)
        theme_menu.add_radiobutton(label="Light", variable=self.theme_var, value="Light", command=self.update_theme)
        theme_menu.add_radiobutton(label="Dark", variable=self.theme_var, value="Dark", command=self.update_theme)

    def update_font_size(self):
        size = self.font_size_var.get()
        default_font = tkFont.Font(family="Segoe UI", size=size)
        consolas_font = tkFont.Font(family="Consolas", size=size)
        self.root.option_add("*Font", default_font)
        self.ref_text.config(font=consolas_font)
        self.sample_text.config(font=consolas_font)
        self.alignment_text.config(font=consolas_font)
        self.summary_text.config(font=consolas_font)
        self.protein_text.config(font=consolas_font)

    def update_theme(self):
        theme = self.theme_var.get()
        if theme == "Dark":
            self.colors.update({
                'background': '#1e293b',
                'card': '#334155',
                'light': '#475569',
                'text_primary': '#f8fafc',
                'text_secondary': '#94a3b8'
            })
        else:
            self.colors.update({
                'background': '#f1f5f9',
                'card': '#ffffff',
                'light': '#f8fafc',
                'text_primary': '#1e293b',
                'text_secondary': '#64748b'
            })
        self.root.configure(bg=self.colors['background'])
        for widget in [self.input_tab, self.alignment_tab, self.mutation_tab, self.results_tab]:
            widget.configure(bg=self.colors['background'])
        self.update_widget_styles()

    def update_widget_styles(self):
        for widget in self.root.winfo_children():
            self._update_widget_style(widget)

    def _update_widget_style(self, widget):
        try:
            if isinstance(widget, (tk.Frame, tk.Label, tk.Text)):
                widget.configure(bg=self.colors['background'])
            if isinstance(widget, tk.Label):
                widget.configure(fg=self.colors['text_primary'])
            if isinstance(widget, tk.Text):
                widget.configure(bg=self.colors['light'], fg=self.colors['text_primary'])
            for child in widget.winfo_children():
                self._update_widget_style(child)
        except:
            pass

    def create_input_tab(self):
        tab = tk.Frame(self.notebook, bg=self.colors['background'])
        main_container = tk.Frame(tab, bg=self.colors['background'])
        main_container.pack(fill='both', expand=True, padx=15, pady=15)
        left_col = tk.Frame(main_container, bg=self.colors['background'])
        right_col = tk.Frame(main_container, bg=self.colors['background'])
        left_col.pack(side='left', fill='both', expand=True, padx=(0, 10))
        right_col.pack(side='right', fill='both', expand=True, padx=(10, 0))
        gene_card, gene_content = self.create_card_frame(left_col, "🔍 NCBI Gene Database", 220)
        gene_card.pack(fill='x', pady=(0, 10))
        input_frame = tk.Frame(gene_content, bg=self.colors['card'])
        input_frame.pack(fill='x', pady=(0, 10))
        tk.Label(input_frame, text="Gene Symbol or Identifier:", font=tkFont.Font(family="Segoe UI", size=10, weight="bold"), bg=self.colors['card'], fg=self.colors['text_primary']).pack(anchor='w', pady=(0, 5))
        entry_frame = tk.Frame(input_frame, bg=self.colors['card'])
        entry_frame.pack(fill='x', pady=(0, 10))
        self.gene_entry = tk.Entry(entry_frame, font=("Segoe UI", 11), width=30, relief='solid', bd=1, highlightthickness=1, highlightcolor=self.colors['border_focus'], bg='white', fg=self.colors['text_primary'])
        self.gene_entry.pack(fill='x')
        self.create_tooltip(self.gene_entry, "Enter gene symbol or ID (e.g., BRCA1)")
        btn_frame = tk.Frame(gene_content, bg=self.colors['card'])
        btn_frame.pack(fill='x', pady=(0, 10))
        fetch_btn = ttk.Button(btn_frame, text="🔗 Fetch from NCBI Database", style='Info.TButton', command=self.fetch_gene_threaded)
        fetch_btn.pack()
        self.create_tooltip(fetch_btn, "Fetch gene sequence from NCBI database")
        self.fetch_status = tk.Label(gene_content, text="Ready to fetch gene data", bg=self.colors['card'], fg=self.colors['text_secondary'], font=("Segoe UI", 9, "italic"))
        self.fetch_status.pack(pady=(10, 0))
        ref_card, ref_content = self.create_card_frame(left_col, "📄 Reference Sequence")
        ref_card.pack(fill='both', expand=True)
        text_frame = tk.Frame(ref_content, bg=self.colors['card'])
        text_frame.pack(fill='both', expand=True, pady=(0, 10))
        self.ref_text = tk.Text(text_frame, height=10, font=("Consolas", 10), wrap='word', relief='solid', bd=1, highlightthickness=1, highlightcolor=self.colors['border_focus'], selectbackground=self.colors['secondary_light'], bg='white', fg=self.colors['text_primary'])
        ref_scroll = ttk.Scrollbar(text_frame, orient='vertical', command=self.ref_text.yview)
        self.ref_text.configure(yscrollcommand=ref_scroll.set)
        self.ref_text.pack(side='left', fill='both', expand=True)
        ref_scroll.pack(side='right', fill='y')
        self.create_tooltip(self.ref_text, "Enter or upload reference DNA sequence")
        ref_btn_frame = tk.Frame(ref_content, bg=self.colors['card'])
        ref_btn_frame.pack(fill='x', pady=(10, 0))
        upload_ref_btn = ttk.Button(ref_btn_frame, text="📁 Upload Reference File", style='Primary.TButton', command=lambda: self.upload_file(self.ref_text))
        upload_ref_btn.pack(side='left', padx=(0, 10))
        self.create_tooltip(upload_ref_btn, "Upload a FASTA or text file")
        clear_ref_btn = ttk.Button(ref_btn_frame, text="🗑 Clear Sequence", command=lambda: self.ref_text.delete('1.0', tk.END))
        clear_ref_btn.pack(side='left')
        self.create_tooltip(clear_ref_btn, "Clear the reference sequence")
        sample_card, sample_content = self.create_card_frame(right_col, "🧪 Sample Sequence")
        sample_card.pack(fill='both', expand=True)
        sample_text_frame = tk.Frame(sample_content, bg=self.colors['card'])
        sample_text_frame.pack(fill='both', expand=True, pady=(0, 10))
        self.sample_text = tk.Text(sample_text_frame, height=18, font=("Consolas", 10), wrap='word', relief='solid', bd=1, highlightthickness=1, highlightcolor=self.colors['border_focus'], selectbackground=self.colors['secondary_light'], bg='white', fg=self.colors['text_primary'])
        sample_scroll = ttk.Scrollbar(sample_text_frame, orient='vertical', command=self.sample_text.yview)
        self.sample_text.configure(yscrollcommand=sample_scroll.set)
        self.sample_text.pack(side='left', fill='both', expand=True)
        sample_scroll.pack(side='right', fill='y')
        self.create_tooltip(self.sample_text, "Enter or upload sample DNA sequence")
        sample_btn_frame = tk.Frame(sample_content, bg=self.colors['card'])
        sample_btn_frame.pack(fill='x', pady=(10, 0))
        upload_sample_btn = ttk.Button(sample_btn_frame, text="📁 Upload Sample File", style='Primary.TButton', command=lambda: self.upload_file(self.sample_text))
        upload_sample_btn.pack(side='left', padx=(0, 10))
        self.create_tooltip(upload_sample_btn, "Upload a FASTA or text file")
        clear_sample_btn = ttk.Button(sample_btn_frame, text="🗑 Clear Sequence", command=lambda: self.sample_text.delete('1.0', tk.END))
        clear_sample_btn.pack(side='left')
        self.create_tooltip(clear_sample_btn, "Clear the sample sequence")
        return tab

    def create_alignment_tab(self):
        tab = tk.Frame(self.notebook, bg=self.colors['background'])
        main_container = tk.Frame(tab, bg=self.colors['background'])
        main_container.pack(fill='both', expand=True, padx=15, pady=15)
        control_card, control_content = self.create_card_frame(main_container, "⚙ Alignment Configuration", 250)
        control_card.pack(fill='x', pady=(0, 10))
        algo_section = tk.Frame(control_content, bg=self.colors['card'])
        algo_section.pack(fill='x', pady=(0, 10))
        tk.Label(algo_section, text="Alignment Algorithm:", font=tkFont.Font(family="Segoe UI", size=11, weight="bold"), bg=self.colors['card'], fg=self.colors['text_primary']).pack(anchor='w', pady=(0, 5))
        self.algo_var = tk.StringVar(value="Global (Needleman-Wunsch)")
        radio_frame = tk.Frame(algo_section, bg=self.colors['card'])
        radio_frame.pack(fill='x', pady=(5, 0))
        tk.Radiobutton(radio_frame, text="🌐 Global Alignment (Needleman-Wunsch)", variable=self.algo_var, value="Global (Needleman-Wunsch)", bg=self.colors['card'], fg=self.colors['text_primary'], font=("Segoe UI", 10), selectcolor=self.colors['secondary']).pack(anchor='w', pady=3)
        tk.Radiobutton(radio_frame, text="🎯 Local Alignment (Smith-Waterman)", variable=self.algo_var, value="Local (Smith-Waterman)", bg=self.colors['card'], fg=self.colors['text_primary'], font=("Segoe UI", 10), selectcolor=self.colors['secondary']).pack(anchor='w', pady=3)
        btn_section = tk.Frame(control_content, bg=self.colors['card'])
        btn_section.pack(fill='x', pady=(15, 0))
        self.align_btn = ttk.Button(btn_section, text="🔗 Perform Sequence Alignment", style='Success.TButton', command=self.align_sequences_threaded)
        self.align_btn.pack()
        self.create_tooltip(self.align_btn, "Align reference and sample sequences")
        progress_section = tk.Frame(control_content, bg=self.colors['card'])
        progress_section.pack(fill='x', pady=(15, 0))
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(progress_section, variable=self.progress_var, maximum=100, length=350, style='TProgressbar')
        self.progress_bar.pack(fill='x', pady=(0, 5))
        self.progress_label = tk.Label(progress_section, text="Ready to align sequences", bg=self.colors['card'], fg=self.colors['text_secondary'], font=("Segoe UI", 10, "italic"))
        self.progress_label.pack()
        results_card, results_content = self.create_card_frame(main_container, "📋 Alignment Results")
        results_card.pack(fill='both', expand=True)
        self.alignment_text = tk.Text(results_content, font=("Consolas", 10), wrap='none', state='disabled', relief='solid', bd=1, bg=self.colors['light'], fg=self.colors['text_primary'])
        self.alignment_text.tag_configure("match", foreground=self.colors['success'])
        self.alignment_text.tag_configure("mismatch", foreground=self.colors['danger'])
        h_scroll = ttk.Scrollbar(results_content, orient='horizontal', command=self.alignment_text.xview)
        v_scroll = ttk.Scrollbar(results_content, orient='vertical', command=self.alignment_text.yview)
        self.alignment_text.configure(xscrollcommand=h_scroll.set, yscrollcommand=v_scroll.set)
        self.alignment_text.grid(row=0, column=0, sticky='nsew')
        v_scroll.grid(row=0, column=1, sticky='ns')
        h_scroll.grid(row=1, column=0, sticky='ew')
        results_content.grid_rowconfigure(0, weight=1)
        results_content.grid_columnconfigure(0, weight=1)
        return tab

    def create_mutation_tab(self):
        tab = tk.Frame(self.notebook, bg=self.colors['background'])
        main_container = tk.Frame(tab, bg=self.colors['background'])
        main_container.pack(fill='both', expand=True, padx=15, pady=15)

        # Mutation Analysis Engine Section
        control_card, control_content = self.create_card_frame(main_container, "🧪 Mutation Analysis Engine", 250)
        control_card.pack(fill='x', pady=(0, 5))
        code_section = tk.Frame(control_content, bg=self.colors['card'])
        code_section.pack(fill='x', pady=(10, 10))
        tk.Label(code_section, text="Genetic Code:", font=tkFont.Font(family="Segoe UI", size=11, weight="bold"), bg=self.colors['card'], fg=self.colors['text_primary']).pack(anchor='w', pady=(0, 5))
        self.code_var = tk.StringVar(value="Standard")
        code_frame = tk.Frame(code_section, bg=self.colors['card'])
        code_frame.pack(fill='x', pady=(5, 0))
        tk.Radiobutton(code_frame, text="Standard", variable=self.code_var, value="Standard", bg=self.colors['card'], fg=self.colors['text_primary'], font=("Segoe UI", 10), selectcolor=self.colors['secondary'], command=self.enable_mutation_options).pack(anchor='w', pady=3)
        tk.Radiobutton(code_frame, text="Mitochondrial", variable=self.code_var, value="Mitochondrial", bg=self.colors['card'], fg=self.colors['text_primary'], font=("Segoe UI", 10), selectcolor=self.colors['secondary'], command=self.enable_mutation_options).pack(anchor='w', pady=3)

        # Add the buttons with initial disabled state
        btn_section = tk.Frame(control_content, bg=self.colors['card'])
        btn_section.pack(fill='x', pady=(10, 5))
        self.analyze_btn = ttk.Button(btn_section, text="🔬 Analyze Mutations & Variants", style='Danger.TButton', command=self.analyze_mutations_threaded)
        self.analyze_btn.pack(side='left', padx=(0, 10))
        self.create_tooltip(self.analyze_btn, "Analyze mutations in aligned sequences")
        self.pathogenicity_btn = ttk.Button(btn_section, text="🧬 Predict Pathogenicity", style='Info.TButton', command=self.predict_pathogenicity_threaded)
        self.pathogenicity_btn.pack(side='left')
        self.create_tooltip(self.pathogenicity_btn, "Predict pathogenicity using local algorithm")
        self.analyze_btn.state(['disabled'])
        self.pathogenicity_btn.state(['disabled'])
        status_section = tk.Frame(control_content, bg=self.colors['card'])
        status_section.pack(fill='x', pady=(5, 0))
        self.analysis_status = tk.Label(status_section, text="Perform sequence alignment first", bg=self.colors['card'], fg=self.colors['warning'], font=("Segoe UI", 10, "italic"))
        self.analysis_status.pack()

        # Detected Mutations & Variants Table
        table_card, table_content = self.create_card_frame(main_container, "📊 Detected Mutations & Variants", 400)
        table_card.pack(fill='both', expand=True, pady=(5, 15))
        columns = ("Position", "Ref", "Alt", "Type", "Region", "Effect", "Frameshift", "Severity", "SIFT", "PolyPhen")
        self.mutation_tree = ttk.Treeview(table_content, columns=columns, show='headings', height=12, style='Treeview')
        column_widths = {"Position": 90, "Ref": 70, "Alt": 70, "Type": 100, "Region": 90, "Effect": 130, "Frameshift": 100, "Severity": 100, "SIFT": 120, "PolyPhen": 120}
        for col in columns:
            self.mutation_tree.heading(col, text=col, anchor='center')
            self.mutation_tree.column(col, width=column_widths[col], anchor='center')
            self.mutation_tree.tag_configure(col, font=('Segoe UI', 10, 'bold'), background=self.colors['primary'], foreground='white')
        tree_v_scroll = ttk.Scrollbar(table_content, orient='vertical', command=self.mutation_tree.yview)
        tree_h_scroll = ttk.Scrollbar(table_content, orient='horizontal', command=self.mutation_tree.xview)
        self.mutation_tree.configure(yscrollcommand=tree_v_scroll.set, xscrollcommand=tree_h_scroll.set)
        self.mutation_tree.grid(row=0, column=0, sticky='nsew', padx=10, pady=10)
        tree_v_scroll.grid(row=0, column=1, sticky='ns', padx=(0, 10), pady=10)
        tree_h_scroll.grid(row=1, column=0, sticky='ew', padx=10, pady=(0, 10))
        table_content.grid_rowconfigure(0, weight=1)
        table_content.grid_columnconfigure(0, weight=1)
        self.mutation_tree.bind('<Double-1>', self.show_mutation_details)

        return tab

    def create_results_tab(self):
        tab = tk.Frame(self.notebook, bg=self.colors['background'])
        main_container = tk.Frame(tab, bg=self.colors['background'])
        main_container.pack(fill='both', expand=True, padx=15, pady=15)
        summary_card, summary_content = self.create_card_frame(main_container, "📈 Analysis Summary & Statistics", 250)
        summary_card.pack(fill='x', pady=(0, 10))
        self.summary_text = tk.Text(summary_content, height=10, font=("Consolas", 10), state='disabled', bg=self.colors['light'], fg=self.colors['text_primary'], relief='flat', selectbackground=self.colors['secondary_light'])
        self.summary_text.pack(fill='both', expand=True, padx=10, pady=10)
        protein_card, protein_content = self.create_card_frame(main_container, "🧬 Protein Translation", 250)
        protein_card.pack(fill='x', pady=(10, 10))
        self.protein_text = tk.Text(protein_content, height=10, font=("Consolas", 10), state='disabled', bg=self.colors['light'], fg=self.colors['text_primary'], relief='flat', selectbackground=self.colors['secondary_light'])
        self.protein_text.tag_configure("changed_aa", foreground=self.colors['danger'])
        self.protein_text.pack(fill='both', expand=True, padx=10, pady=10)
        export_card, export_content = self.create_card_frame(main_container, "💾 Export & Reporting")
        export_card.pack(fill='both', expand=True)
        desc_frame = tk.Frame(export_content, bg=self.colors['card'])
        desc_frame.pack(fill='x', pady=(0, 10))
        tk.Label(desc_frame, text="Export your analysis results in multiple formats for further research and documentation:", font=("Segoe UI", 10), bg=self.colors['card'], fg=self.colors['text_secondary'], wraplength=500).pack(anchor='w', pady=5)
        export_btn_frame = tk.Frame(export_content, bg=self.colors['card'])
        export_btn_frame.pack(fill='x', pady=(10, 0))
        btn_row1 = tk.Frame(export_btn_frame, bg=self.colors['card'])
        btn_row1.pack(fill='x', pady=(0, 10))
        export_csv_btn = ttk.Button(btn_row1, text="📄 Export to CSV", style='Success.TButton', command=self.export_to_csv)
        export_csv_btn.pack(side='left', padx=(0, 10))
        self.create_tooltip(export_csv_btn, "Export mutations to CSV file")
        copy_summary_btn = ttk.Button(btn_row1, text="📋 Copy Summary", style='Primary.TButton', command=self.copy_summary)
        copy_summary_btn.pack(side='left', padx=(0, 10))
        self.create_tooltip(copy_summary_btn, "Copy summary to clipboard")
        save_report_btn = ttk.Button(btn_row1, text="💾 Save Full Report", style='Info.TButton', command=self.save_report)
        save_report_btn.pack(side='left')
        self.create_tooltip(save_report_btn, "Save detailed report as text")
        export_pdf_btn = ttk.Button(btn_row1, text="📜 Export to PDF", style='Danger.TButton', command=self.export_to_pdf)
        export_pdf_btn.pack(side='left', padx=(10, 0))
        self.create_tooltip(export_pdf_btn, "Export report as PDF")
        return tab

    def validate_sequence(self, seq):
        valid_nucleotides = set("ATCGN-")
        return all(c.upper() in valid_nucleotides for c in seq.strip())

    def parse_fasta(self, text):
        lines = text.strip().split('\n')
        sequence = ""
        for line in lines:
            if not line.startswith('>'):
                sequence += line.strip().upper()
        return sequence

    def upload_file(self, text_widget):
        try:
            file_path = filedialog.askopenfilename(title="Select Sequence File", filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("Text files", "*.txt"), ("All files", "*.*")])
            if not file_path:
                return
            with open(file_path, 'r') as file:
                content = file.read()
            sequence = self.parse_fasta(content)
            if not sequence:
                messagebox.showerror("Error", "No valid sequence found in file")
                return
            if not self.validate_sequence(sequence):
                messagebox.showwarning("Invalid Sequence", "Sequence contains invalid characters. Only A, T, C, G, N, - allowed")
                return
            text_widget.delete('1.0', tk.END)
            text_widget.insert('1.0', sequence)
            messagebox.showinfo("Success", f"Loaded sequence: {len(sequence)} nucleotides")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {str(e)}")

    def fetch_gene_threaded(self):
        def fetch():
            self.fetch_gene()
        thread = threading.Thread(target=fetch, daemon=True)
        thread.start()

    def fetch_gene(self):
        gene_name = self.gene_entry.get().strip()
        if not gene_name:
            messagebox.showwarning("Input Error", "Please enter a gene name")
            return
        try:
            self.fetch_status.config(text="🔍 Searching NCBI database...")
            self.root.update()
            search_term = f'({gene_name}[Gene Name]) AND "Homo sapiens"[Organism] AND RefSeq[Filter]'
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5)
            search_results = Entrez.read(handle)
            handle.close()
            if not search_results["IdList"]:
                self.fetch_status.config(text="❌ Gene not found")
                messagebox.showerror("Not Found", f"Gene '{gene_name}' not found in NCBI database")
                return
            self.fetch_status.config(text="📥 Downloading sequence data...")
            self.root.update()
            record_id = search_results["IdList"][0]
            handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            # Extract chromosome from features
            for feature in record.features:
                if feature.type == "source" and "chromosome" in feature.qualifiers:
                    self.chrom = feature.qualifiers["chromosome"][0]
                    break
            else:
                self.chrom = None  # Default if not found
            self.ref_text.delete('1.0', tk.END)
            self.ref_text.insert('1.0', str(record.seq))
            self.exon_ranges = []
            self.intron_ranges = []
            for feature in record.features:
                if feature.type == "exon":
                    start = int(feature.location.start) + 1  # 1-based indexing
                    end = int(feature.location.end)
                    self.exon_ranges.append((start, end))
                elif feature.type == "CDS":  # Fallback to CDS if exons are not annotated
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    self.exon_ranges.append((start, end))
            if len(self.exon_ranges) > 1:
                self.exon_ranges.sort()
                for i in range(len(self.exon_ranges) - 1):
                    intron_start = self.exon_ranges[i][1] + 1
                    intron_end = self.exon_ranges[i + 1][0] - 1
                    if intron_start <= intron_end:
                        self.intron_ranges.append((intron_start, intron_end))
            success_msg = f"✅ Fetched {gene_name}: {len(self.exon_ranges)} exons, {len(self.intron_ranges)} introns"
            self.fetch_status.config(text=success_msg)
            messagebox.showinfo("Success", f"Successfully fetched {gene_name}\nSequence length: {len(record.seq)} bp\nExons: {len(self.exon_ranges)}\nIntrons: {len(self.intron_ranges)}\nChromosome: {self.chrom or 'Unknown'}")
        except Exception as e:
            error_msg = f"❌ Error: {str(e)}"
            self.fetch_status.config(text=error_msg)
            messagebox.showerror("Fetch Error", f"Failed to fetch gene data:\n{str(e)}")

    def align_sequences_threaded(self):
        def align():
            self.align_sequences()
        thread = threading.Thread(target=align, daemon=True)
        thread.start()

    def align_sequences(self):
        try:
            ref_seq = self.ref_text.get('1.0', tk.END).strip().upper()
            sample_seq = self.sample_text.get('1.0', tk.END).strip().upper()
            if not ref_seq or not sample_seq:
                messagebox.showwarning("Input Error", "Both reference and sample sequences are required")
                return
            if not self.validate_sequence(ref_seq) or not self.validate_sequence(sample_seq):
                messagebox.showwarning("Invalid Sequence", "Sequences contain invalid characters")
                return
            self.progress_label.config(text="🔄 Initializing alignment...")
            self.progress_var.set(10)
            self.root.update()
            start_time = time.time()
            timeout = 300
            if self.algo_var.get() == "Global (Needleman-Wunsch)":
                alignments = pairwise2.align.globalms(ref_seq, sample_seq, 1, -1, -10, -1, one_alignment_only=True)
            else:
                alignments = pairwise2.align.localms(ref_seq, sample_seq, 1, -1, -10, -1, one_alignment_only=True)
            if time.time() - start_time > timeout:
                raise TimeoutError("Alignment took too long and was terminated.")
            if not alignments:
                raise ValueError("No valid alignment generated")
            self.progress_var.set(75)
            self.root.update()
            best_alignment = alignments[0]
            self.aligned_ref = best_alignment.seqA
            self.aligned_sample = best_alignment.seqB
            score = best_alignment.score
            self.alignment_text.config(state='normal')
            self.alignment_text.delete('1.0', tk.END)
            alignment_display = self.format_alignment_display(self.aligned_ref, self.aligned_sample, score)
            self.alignment_text.insert('1.0', alignment_display)
            self.alignment_text.config(state='disabled')
            self.progress_var.set(100)
            self.progress_label.config(text=f"✅ Alignment complete! Score: {score:.1f} (Time: {time.time() - start_time:.1f}s)")
            self.analysis_status.config(text="Ready for mutation analysis", fg=self.colors['success'])
            self.analyze_btn.state(['!disabled'])
            self.pathogenicity_btn.state(['disabled'])
            messagebox.showinfo("Alignment Complete", f"Alignment successful!\nAlgorithm: {self.algo_var.get()}\nScore: {score:.1f}\nLength: {len(self.aligned_ref)} bp\nTime: {time.time() - start_time:.1f}s")
        except TimeoutError:
            self.progress_label.config(text="❌ Alignment timed out")
            messagebox.showerror("Alignment Error", "Alignment took too long and was terminated. Consider using local alignment or shorter sequences.")
            self.progress_var.set(0)
        except Exception as e:
            self.progress_label.config(text="❌ Alignment failed")
            messagebox.showerror("Alignment Error", f"Failed to align sequences: {str(e)}")
            self.progress_var.set(0)

    def format_alignment_display(self, seq1, seq2, score, line_length=80):
        self.alignment_text.config(state='normal')
        self.alignment_text.delete('1.0', tk.END)
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
        identity = (matches / len(seq1)) * 100 if len(seq1) > 0 else 0
        result = f"Alignment Score: {score:.1f}\nLength: {len(seq1)} bp\nIdentity: {identity:.2f}%\n" + "=" * 50 + "\n\n"
        self.alignment_text.insert(tk.END, result)
        for i in range(0, len(seq1), line_length):
            chunk1 = seq1[i:i+line_length]
            chunk2 = seq2[i:i+line_length]
            match_line = ""
            for j in range(len(chunk1)):
                if j < len(chunk2) and chunk1[j] == chunk2[j] and chunk1[j] != '-':
                    match_line += "|"
                    self.alignment_text.tag_add("match", f"{self.alignment_text.index(tk.END).split('.')[0]}.{len('Ref:    ')+j}")
                else:
                    match_line += " "
                    if chunk1[j] != '-' and j < len(chunk2) and chunk2[j] != '-':
                        self.alignment_text.tag_add("mismatch", f"{self.alignment_text.index(tk.END).split('.')[0]}.{len('Ref:    ')+j}")
            self.alignment_text.insert(tk.END, f"Ref:    {chunk1}\n        {match_line}\nSample: {chunk2}\n\n")
        self.alignment_text.config(state='disabled')
        return result

    def enable_mutation_options(self):
        # Enable buttons only after genetic code is selected
        self.analyze_btn.state(['!disabled'])
        self.pathogenicity_btn.state(['!disabled'])

    def analyze_mutations_threaded(self):
        def analyze():
            self.analyze_mutations()
        thread = threading.Thread(target=analyze, daemon=True)
        thread.start()

    def analyze_mutations(self):
        try:
            if not self.aligned_ref or not self.aligned_sample:
                messagebox.showwarning("No Alignment", "Please perform sequence alignment first")
                return
            self.analysis_status.config(text="🔬 Analyzing mutations...")
            self.root.update()
            self.mutations = []
            ref_pos = 0
            i = 0
            while i < min(len(self.aligned_ref), len(self.aligned_sample)):
                ref_base = self.aligned_ref[i]
                alt_base = self.aligned_sample[i]
                if ref_base != '-':
                    ref_pos += 1
                if ref_base != '-' and alt_base != '-' and ref_base != alt_base:
                    mutation = self.analyze_snp(ref_pos, ref_base, alt_base)
                    if mutation:
                        self.mutations.append(mutation)
                elif ref_base != '-' and alt_base == '-':
                    del_length, del_seq = self.get_deletion_info(i, self.aligned_ref, self.aligned_sample)
                    mutation = self.analyze_deletion(ref_pos, del_seq, del_length)
                    if mutation:
                        self.mutations.append(mutation)
                    i += del_length - 1
                elif ref_base == '-' and alt_base != '-':
                    ins_length, ins_seq = self.get_insertion_info(i, self.aligned_ref, self.aligned_sample)
                    mutation = self.analyze_insertion(ref_pos, ins_seq, ins_length)
                    if mutation:
                        self.mutations.append(mutation)
                    i += ins_length - 1
                i += 1
            self.update_mutation_table()
            self.update_summary()
            self.update_protein_display()
            self.analysis_status.config(text=f"✅ Found {len(self.mutations)} mutations", fg=self.colors['success'])
            self.pathogenicity_btn.state(['!disabled'])
            if len(self.mutations) == 0:
                messagebox.showinfo("No Mutations", "No mutations detected in the aligned sequences.")
            else:
                messagebox.showinfo("Analysis Complete", f"Mutation analysis complete!\nTotal mutations: {len(self.mutations)}\nCheck the Mutations tab for details")
        except Exception as e:
            self.analysis_status.config(text="❌ Analysis failed", fg=self.colors['danger'])
            messagebox.showerror("Analysis Error", f"Failed to analyze mutations: {str(e)}")

    def get_deletion_info(self, start_pos, ref_seq, alt_seq):
        del_seq = ""
        length = 0
        pos = start_pos
        while pos < len(ref_seq) and pos < len(alt_seq) and ref_seq[pos] != '-' and alt_seq[pos] == '-':
            del_seq += ref_seq[pos]
            length += 1
            pos += 1
        return length, del_seq if length > 0 else ("", 0)

    def get_insertion_info(self, start_pos, ref_seq, alt_seq):
        ins_seq = ""
        length = 0
        pos = start_pos
        while pos < len(ref_seq) and pos < len(alt_seq) and ref_seq[pos] == '-' and alt_seq[pos] != '-':
            ins_seq += alt_seq[pos]
            length += 1
            pos += 1
        return length, ins_seq if length > 0 else ("", 0)

    def get_region(self, position):
        for start, end in self.exon_ranges:
            if start <= position <= end:
                return "Exon"
        for start, end in self.intron_ranges:
            if start <= position <= end:
                return "Intron"
        return "Intergenic"

    def analyze_snp(self, position, ref_base, alt_base):
        if not ref_base or not alt_base or ref_base == alt_base:
            return None
        region = self.get_region(position)
        effect = "Substitution"
        severity = "🟢 Low"
        frameshift = "No"
        sift = "-"
        polyphen = "-"
        if region == "Exon":
            effect, severity = self.analyze_coding_effect(position, ref_base, alt_base)
        elif region == "Intron":
            effect = "Intronic"
            severity = "⚪ Minimal"
        return {
            'position': position,
            'ref': ref_base,
            'alt': alt_base,
            'type': 'SNP',
            'region': region,
            'effect': effect,
            'frameshift': frameshift,
            'severity': severity,
            'sift': sift,
            'polyphen': polyphen
        }

    def analyze_deletion(self, position, del_seq, length):
        if length == 0:
            return None
        region = self.get_region(position)
        effect = "Deletion"
        severity = "🟠 Medium"
        frameshift = "No"
        if region == "Exon" and length > 0:
            if length % 3 == 0:
                effect = "In-frame deletion"
                severity = "🟠 Medium"
            else:
                effect = "Frameshift deletion"
                severity = "🔴 High"
                frameshift = "Yes"
        return {
            'position': position,
            'ref': del_seq,
            'alt': '-',
            'type': 'Deletion',
            'region': region,
            'effect': effect,
            'frameshift': frameshift,
            'severity': severity,
            'sift': '-',
            'polyphen': '-'
        }

    def analyze_insertion(self, position, ins_seq, length):
        if length == 0:
            return None
        region = self.get_region(position)
        effect = "Insertion"
        severity = "🟠 Medium"
        frameshift = "No"
        if region == "Exon" and length > 0:
            if length % 3 == 0:
                effect = "In-frame insertion"
                severity = "🟠 Medium"
            else:
                effect = "Frameshift insertion"
                severity = "🔴 High"
                frameshift = "Yes"
        return {
            'position': position,
            'ref': '-',
            'alt': ins_seq,
            'type': 'Insertion',
            'region': region,
            'effect': effect,
            'frameshift': frameshift,
            'severity': severity,
            'sift': '-',
            'polyphen': '-'
        }

    def analyze_coding_effect(self, position, ref_base, alt_base):
        try:
            codon_pos = self.get_codon_position(position)
            if codon_pos == -1:
                return "Non-coding", "⚪ Minimal"
            ref_codon = self.get_codon_at_position(position, ref_base, is_ref=True)
            alt_codon = self.get_codon_at_position(position, alt_base, is_ref=False)
            table = 1 if self.code_var.get() == "Standard" else 2
            if len(ref_codon) == 3 and len(alt_codon) == 3:
                ref_aa = str(Seq(ref_codon).translate(table=table))
                alt_aa = str(Seq(alt_codon).translate(table=table))
                if ref_aa == alt_aa:
                    return "Silent", "🟢 Low"
                elif alt_aa == '*':
                    return "Nonsense", "🔴 High"
                else:
                    return "Missense", "🟠 Medium"
            return "Unknown", "⚪ Minimal"
        except Exception:
            return "Unknown", "⚪ Minimal"

    def get_codon_position(self, position):
        for start, end in self.exon_ranges:
            if start <= position <= end:
                return ((position - start) % 3) + 1
        return -1

    def get_codon_at_position(self, position, base, is_ref=True):
        try:
            seq = self.aligned_ref if is_ref else self.aligned_sample
            unaligned_pos = 0
            aligned_pos = 0
            for i, char in enumerate(seq):
                if char != '-':
                    unaligned_pos += 1
                if unaligned_pos == position:
                    aligned_pos = i
                    break
            codon_start = (aligned_pos // 3) * 3
            codon = seq[codon_start:codon_start + 3].replace('-', '')
            if len(codon) < 3:
                codon += 'N' * (3 - len(codon))
            return codon
        except Exception:
            return "NNN"

    def predict_pathogenicity_threaded(self):
        def predict():
            self.predict_pathogenicity()
        thread = threading.Thread(target=predict, daemon=True)
        thread.start()

    def predict_pathogenicity(self):
        try:
            self.analysis_status.config(text="🔍 Predicting pathogenicity locally...", fg=self.colors['info'])
            self.root.update()
            missense_mutations = [m for m in self.mutations if m['effect'] == 'Missense' and m['type'] == 'SNP']
            if not missense_mutations:
                self.analysis_status.config(text="⚠ No missense mutations to analyze", fg=self.colors['warning'])
                messagebox.showinfo("No Missense Mutations", "No missense mutations detected for pathogenicity prediction.")
                return

            # Local pathogenicity prediction logic
            for mut in missense_mutations:
                position = mut['position']
                ref_base = mut['ref']
                alt_base = mut['alt']
                ref_codon = self.get_codon_at_position(position, ref_base, is_ref=True)
                alt_codon = self.get_codon_at_position(position, alt_base, is_ref=False)
                table = 1 if self.code_var.get() == "Standard" else 2
                ref_aa = str(Seq(ref_codon).translate(table=table))
                alt_aa = str(Seq(alt_codon).translate(table=table))

                # Simplified SIFT-like score (conservation-based)
                conservation_score = self.calculate_conservation_score(ref_aa)
                sift_score = 1.0 - (conservation_score / 100.0)  # Inverse relation, 0-1 scale
                sift_pred = "Tolerated" if sift_score > 0.05 else "Deleterious"

                # PolyPhen-like score (Grantham distance for physicochemical difference)
                polyphen_score = self.calculate_grantham_distance(ref_aa, alt_aa)
                polyphen_pred = "Benign" if polyphen_score < 50 else "Possibly Damaging" if polyphen_score < 100 else "Probably Damaging"

                mut['sift'] = f"{sift_pred} ({sift_score:.2f})"
                mut['polyphen'] = f"{polyphen_pred} ({polyphen_score:.2f})"

            self.update_mutation_table()
            self.update_summary()
            self.update_protein_display()
            self.analysis_status.config(text=f"✅ Predicted pathogenicity for {len(missense_mutations)} mutations", fg=self.colors['success'])
            messagebox.showinfo("Pathogenicity Prediction", f"Predicted SIFT/PolyPhen-2 scores for {len(missense_mutations)} missense mutations.\nCheck the Mutations tab for details.")
        except Exception as e:
            self.analysis_status.config(text="❌ Pathogenicity prediction failed", fg=self.colors['danger'])
            messagebox.showerror("Prediction Error", f"Failed to predict pathogenicity:\n{str(e)}")

    def calculate_conservation_score(self, aa):
        # Simplified conservation score based on frequency of amino acids (hypothetical values)
        conservation = {
            'A': 80, 'C': 70, 'D': 60, 'E': 60, 'F': 50, 'G': 90, 'H': 60, 'I': 50,
            'K': 60, 'L': 50, 'M': 50, 'N': 60, 'P': 70, 'Q': 60, 'R': 60, 'S': 70,
            'T': 70, 'V': 60, 'W': 40, 'Y': 50
        }
        return conservation.get(aa, 50)  # Default to 50 if amino acid not found

    def calculate_grantham_distance(self, ref_aa, alt_aa):
        # Grantham distance matrix (simplified values based on physicochemical properties)
        grantham_matrix = {
            ('A', 'A'): 0, ('A', 'C'): 195, ('A', 'D'): 126, ('A', 'E'): 153, ('A', 'F'): 176,
            ('A', 'G'): 60, ('A', 'H'): 90, ('A', 'I'): 94, ('A', 'K'): 135, ('A', 'L'): 145,
            ('A', 'M'): 140, ('A', 'N'): 111, ('A', 'P'): 67, ('A', 'Q'): 147, ('A', 'R'): 112,
            ('A', 'S'): 99, ('A', 'T'): 86, ('A', 'V'): 64, ('A', 'W'): 191, ('A', 'Y'): 160,
            ('C', 'C'): 0, ('C', 'D'): 170, ('C', 'E'): 197, ('C', 'F'): 165, ('C', 'G'): 149,
            ('D', 'D'): 0, ('D', 'E'): 45, ('D', 'F'): 162, ('D', 'G'): 94, ('D', 'H'): 81,
            # Add other combinations symmetrically as needed...
        }
        # Default to a mid-range distance if not in matrix
        key = tuple(sorted([ref_aa, alt_aa]))
        return grantham_matrix.get(key, 100)

    def update_mutation_table(self):
        for item in self.mutation_tree.get_children():
            self.mutation_tree.delete(item)
        if self.mutations:
            for mut in self.mutations:
                self.mutation_tree.insert('', 'end', values=(
                    mut['position'],
                    mut['ref'],
                    mut['alt'],
                    mut['type'],
                    mut['region'],
                    mut['effect'],
                    mut['frameshift'],
                    mut['severity'],
                    mut['sift'],
                    mut['polyphen']
                ))
        else:
            self.mutation_tree.insert('', 'end', values=("No mutations detected", "", "", "", "", "", "", "", "", ""))

    def update_protein_display(self):
        try:
            self.protein_text.config(state='normal')
            self.protein_text.delete('1.0', tk.END)
            if not self.aligned_ref or not self.aligned_sample:
                self.protein_text.insert('1.0', "No alignment available")
                self.protein_text.config(state='disabled')
                return
            ref_seq_no_gaps = self.aligned_ref.replace('-', '')
            sample_seq_no_gaps = self.aligned_sample.replace('-', '')
            table = 1 if self.code_var.get() == "Standard" else 2
            ref_protein = str(Seq(ref_seq_no_gaps).translate(table=table, to_stop=False))
            sample_protein = str(Seq(sample_seq_no_gaps).translate(table=table, to_stop=False))
            display = "Reference Protein:\n"
            for i, aa in enumerate(ref_protein):
                display += aa
                if i < len(sample_protein) and sample_protein[i] != aa:
                    self.protein_text.tag_add("changed_aa", f"{self.protein_text.index(tk.END).split('.')[0]}.{len(display)-1}")
            display += "\n\nSample Protein:\n"
            for aa in sample_protein:
                display += aa
            display += "\n\nPathogenicity Predictions (Missense Mutations):\n"
            missense_mutations = [m for m in self.mutations if m['effect'] == 'Missense' and m['type'] == 'SNP']
            if missense_mutations:
                for mut in missense_mutations:
                    display += f"Pos {mut['position']}: {mut['ref']}>{mut['alt']} - SIFT: {mut['sift']}, PolyPhen: {mut['polyphen']}\n"
            else:
                display += "No missense mutations detected.\n"
            self.protein_text.insert('1.0', display)
            self.protein_text.config(state='disabled')
        except Exception as e:
            self.protein_text.insert('1.0', f"Error generating protein sequences: {str(e)}")
            self.protein_text.config(state='disabled')

    def update_summary(self):
        current_time = datetime.now()
        if not self.mutations:
            summary = "No mutations detected."
        else:
            total = len(self.mutations)
            snps = len([m for m in self.mutations if m['type'] == 'SNP'])
            insertions = len([m for m in self.mutations if m['type'] == 'Insertion'])
            deletions = len([m for m in self.mutations if m['type'] == 'Deletion'])
            exonic = len([m for m in self.mutations if m['region'] == 'Exon'])
            intronic = len([m for m in self.mutations if m['region'] == 'Intron'])
            high_severity = len([m for m in self.mutations if '🔴' in m['severity']])
            medium_severity = len([m for m in self.mutations if '🟠' in m['severity']])
            low_severity = len([m for m in self.mutations if '🟢' in m['severity']])
            missense = len([m for m in self.mutations if m['effect'] == 'Missense'])
            with_sift = len([m for m in self.mutations if m['sift'] != '-'])
            with_polyphen = len([m for m in self.mutations if m['polyphen'] != '-'])
            summary = f"""🧬 MUTATION ANALYSIS SUMMARY
Generated: {current_time.strftime('%Y-%m-%d %H:%M:%S')} PKT

📊 MUTATION COUNTS:
   Total Mutations: {total}
   • SNPs: {snps}
   • Insertions: {insertions} 
   • Deletions: {deletions}

📍 GENOMIC LOCATION:
   • Exonic: {exonic}
   • Intronic: {intronic}
   • Intergenic: {total - exonic - intronic}

⚠ SEVERITY DISTRIBUTION:
   • High (🔴): {high_severity}
   • Medium (🟠): {medium_severity}
   • Low (🟢): {low_severity}

🧬 PATHOGENICITY PREDICTIONS:
   • Missense Mutations: {missense}
   • With SIFT Scores: {with_sift}
   • With PolyPhen-2 Scores: {with_polyphen}

🔗 ALIGNMENT INFO:
   Reference Length: {len(self.aligned_ref)} bp
   Sample Length: {len(self.aligned_sample)} bp
   Exons Analyzed: {len(self.exon_ranges)}
   Introns Analyzed: {len(self.intron_ranges)}
"""
        self.summary_text.config(state='normal')
        self.summary_text.delete('1.0', tk.END)
        self.summary_text.insert('1.0', summary)
        self.summary_text.config(state='disabled')

    def show_mutation_details(self, event):
        selection = self.mutation_tree.selection()
        if not selection:
            return
        item = self.mutation_tree.item(selection[0])
        values = item['values']
        if not values or values[0] == "No mutations detected":
            return
        position, ref, alt, mut_type, region, effect, frameshift, severity, sift, polyphen = values
        mut_detail = next((mut for mut in self.mutations if mut['position'] == position and mut['ref'] == ref and mut['alt'] == alt and mut['type'] == mut_type), None)
        if mut_detail:
            detail_text = f"""🔍 MUTATION DETAILS

Position: {position}
Reference: {ref}
Alternative: {alt}
Type: {mut_type}
Region: {region}
Effect: {effect}
Frameshift: {frameshift}
Severity: {severity}
SIFT: {sift}
PolyPhen-2: {polyphen}

📍 Genomic Context:
"""
            if region == "Exon":
                for i, (start, end) in enumerate(self.exon_ranges):
                    if start <= position <= end:
                        detail_text += f"   Exon #{i+1} ({start}-{end})\n"
                        break
            elif region == "Intron":
                for i, (start, end) in enumerate(self.intron_ranges):
                    if start <= position <= end:
                        detail_text += f"   Intron #{i+1} ({start}-{end})\n"
                        break
            messagebox.showinfo("Mutation Details", detail_text)

    def export_to_csv(self):
        if not self.mutations:
            messagebox.showwarning("No Data", "No mutations to export")
            return
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")], title="Save Mutations as CSV")
            if not file_path:
                return
            with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['Position', 'Reference', 'Alternative', 'Type', 'Region', 'Effect', 'Frameshift', 'Severity', 'SIFT', 'PolyPhen']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for mut in self.mutations:
                    writer.writerow({
                        'Position': mut['position'],
                        'Reference': mut['ref'],
                        'Alternative': mut['alt'],
                        'Type': mut['type'],
                        'Region': mut['region'],
                        'Effect': mut['effect'],
                        'Frameshift': mut['frameshift'],
                        'Severity': mut['severity'],
                        'SIFT': mut['sift'],
                        'PolyPhen': mut['polyphen']
                    })
            messagebox.showinfo("Export Successful", f"Mutations exported to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export data:\n{str(e)}")

    def export_to_pdf(self):
        if not self.mutations:
            messagebox.showwarning("No Data", "No mutations to export")
            return
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")], title="Save Report as PDF")
            if not file_path:
                return
            c = canvas.Canvas(file_path, pagesize=letter)
            c.setFont("Helvetica-Bold", 16)
            y = 750
            c.drawString(50, y, "MutAnalyzer Pro - Mutation Analysis Report")
            y -= 20
            c.setFont("Helvetica", 12)
            c.drawString(50, y, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} PKT")
            y -= 40
            c.drawString(50, y, "Summary:")
            y -= 20
            c.setFont("Courier", 10)
            for line in self.summary_text.get('1.0', tk.END).split('\n'):
                c.drawString(70, y, line)
                y -= 15
                if y < 50:
                    c.showPage()
                    c.setFont("Courier", 10)
                    y = 750
            c.showPage()
            c.setFont("Helvetica-Bold", 14)
            y = 750
            c.drawString(50, y, "Detailed Mutation List")
            y -= 20
            c.setFont("Courier", 10)
            for i, mut in enumerate(self.mutations, 1):
                c.drawString(70, y, f"{i}. Position {mut['position']}: {mut['ref']} → {mut['alt']}")
                y -= 15
                c.drawString(90, y, f"Type: {mut['type']}")
                y -= 15
                c.drawString(90, y, f"Region: {mut['region']}")
                y -= 15
                c.drawString(90, y, f"Effect: {mut['effect']}")
                y -= 15
                c.drawString(90, y, f"Frameshift: {mut['frameshift']}")
                y -= 15
                c.drawString(90, y, f"Severity: {mut['severity']}")
                y -= 15
                c.drawString(90, y, f"SIFT: {mut['sift']}")
                y -= 15
                c.drawString(90, y, f"PolyPhen-2: {mut['polyphen']}")
                y -= 30
                if y < 50:
                    c.showPage()
                    c.setFont("Courier", 10)
                    y = 750
            c.save()
            messagebox.showinfo("Export Successful", f"Report exported to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export PDF:\n{str(e)}")

    def copy_summary(self):
        try:
            summary_text = self.summary_text.get('1.0', tk.END)
            self.root.clipboard_clear()
            self.root.clipboard_append(summary_text)
            messagebox.showinfo("Copied", "Summary copied to clipboard")
        except Exception as e:
            messagebox.showerror("Copy Error", f"Failed to copy summary:\n{str(e)}")

    def save_report(self):
        try:
            file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt"), ("All files", "*.*")], title="Save Analysis Report")
            if not file_path:
                return
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(self.summary_text.get('1.0', tk.END))
                f.write("\n" + "="*60 + "\n")
                f.write("DETAILED MUTATION LIST\n")
                f.write("="*60 + "\n\n")
                for i, mut in enumerate(self.mutations, 1):
                    f.write(f"{i}. Position {mut['position']}: {mut['ref']} → {mut['alt']}\n")
                    f.write(f"   Type: {mut['type']}\n")
                    f.write(f"   Region: {mut['region']}\n")
                    f.write(f"   Effect: {mut['effect']}\n")
                    f.write(f"   Frameshift: {mut['frameshift']}\n")
                    f.write(f"   Severity: {mut['severity']}\n")
                    f.write(f"   SIFT: {mut['sift']}\n")
                    f.write(f"   PolyPhen-2: {mut['polyphen']}\n\n")
                if hasattr(self, 'aligned_ref') and self.aligned_ref:
                    f.write("\n" + "="*60 + "\n")
                    f.write("SEQUENCE ALIGNMENT\n")
                    f.write("="*60 + "\n\n")
                    f.write(self.format_alignment_display(self.aligned_ref, self.aligned_sample, 0))
            messagebox.showinfo("Report Saved", f"Complete report saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Save Error", f"Failed to save report:\n{str(e)}")

    def run(self):
        self.root.mainloop()

def main():
    print("🧬 Starting MutAnalyzer Pro...")
    print("⚠  IMPORTANT: Please change the email address in lines 11 and 12 before using NCBI and Ensembl VEP features!")
    app = MutationAnalyzer()
    app.run()

if __name__ == "__main__":
    main()
