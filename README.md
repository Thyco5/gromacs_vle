\documentclass{article}
\usepackage{geometry}
\usepackage[svgnames]{xcolor}
\usepackage{titlesec}
\usepackage{listings}
\usepackage{graphicx}
\usepackage{tcolorbox}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage{amsmath} % For math symbols like \rightarrow, \times, \rho, subscripts

% Adjust geometry if needed for more text
\geometry{a4paper, margin=1in}

% --- Your Custom Formatting (from original) ---
\titleformat{\section}{\normalfont\Large\bfseries}{\thesection}{1em}{}
\titlespacing*{\section}{0pt}{2em}{0.5em}
\titleformat{\subsection}{\normalfont\large\bfseries}{\thesubsection}{1em}{}
\titlespacing*{\subsection}{0pt}{1.5em}{0.4em} % Adjusted spacing slightly
\titleformat{\subsubsection}{\normalfont\normalsize\bfseries}{\thesubsubsection}{1em}{}
\titlespacing*{\subsubsection}{0pt}{1em}{0.3em} % Adjusted spacing slightly

% Custom showybox environment
\newtcolorbox{showybox}[2]{
    colframe=#1,
    colback=#1!10,
    title=#2,
    fonttitle=\bfseries\large, % Made title font a bit larger
    arc=3mm, % Slightly more rounded arc
    boxrule=1pt,
    top=6pt, % Increased padding
    bottom=6pt,
    left=6pt,
    right=6pt,
    enlarge top by=0.1cm, % Ensure space for title if it wraps
    breakable, % Allows box to break across pages if content is long
    pad at break=2mm, % Padding at page break
}

% Code listing setup
\lstset{
    basicstyle=\ttfamily\small,
    breaklines=true,
    frame=tb, % Top and bottom frame lines
    framesep=5pt,
    framextopmargin=3pt,
    framexbottommargin=3pt,
    backgroundcolor=\color{LightGray!20},
    captionpos=b,
    aboveskip=1.5em, % Increased spacing
    belowskip=1.5em,
    keywordstyle=\color{blue}\bfseries,
    commentstyle=\color{ForestGreen}\itshape,
    stringstyle=\color{BrickRed},
    numbers=left,
    numberstyle=\tiny\color{DarkGray},
    numbersep=5pt,
    tabsize=2,
    showstringspaces=false,
    extendedchars=true,
    morekeywords={antechamber, parmchk2, packmol, tleap, gmx, grompp, mdrun, editconf, InterMol, python, nohup, loadmol2, loadamberparams, source, set, saveamberparm, quit, sqm, tolerance, filetype, output, structure, number, inside, cube, end, integrator, emtol, emstep, nsteps, cutoff-scheme, nstlist, rlist, coulombtype, rcoulomb, rvdw, pbc, steep, Verlet, PME, xyz, dt, nstxout-compressed, nstvout, nstenergy, nstlog, pme_order, fourierspacing, tcoupl, tc-grps, tau_t, ref_t, gen_vel, gen_temp, gen_seed, constraints, constraint_algorithm, lincs_iter, lincs_order, md, V-rescale, h-bonds, lincs, pcoupltype, tau_p, ref_p, compressibility, Berendsen, isotropic, nstfout, fourier_nx, fourier_ny, fourier_nz, plot, title, xlabel, ylabel, grid, yrange, terminal, pngcairo, unset}, % Added more GROMACS MDP and Gnuplot keywords
    emph={MOH, PRH, GAFF2, AM1-BCC, PDB, MOL2, FRCMOD, PRMTOP, INPCRD, TOP, GRO, MDP, XTC, EDR, TPR, CPT, LINCS, PME, Berendsen, V-rescale, Parrinello-Rahman, NVT, NPT, VLE, UFF, MMFF94, HF/6-31G*, Angstroms, CRYST1, POL, OL, HOP, ParmEd}, % Emphasize key terms
    emphstyle=\color{DarkOrchid}\bfseries,
    literate={->}{$\rightarrow$}{2} {<-}{$\leftarrow$}{2} {>=}{$\geq$}{2} {<=}{$\leq$}{2} % Handle arrows and comparison
}

% Document metadata
\title{Comprehensive Guide: Simulating Vapor-Liquid Equilibrium of a Binary Mixture using AmberTools \& GROMACS}
\author{Mohammad Torikh}
\date{June 21, 2025} % Or \date{\today} for current date

\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,
    urlcolor=teal,
    pdftitle={Comprehensive Guide: Simulating Vapor-Liquid Equilibrium of a Binary Mixture using AmberTools \& GROMACS},
    pdfauthor={Mohammad Torikh},
    pdfsubject={VLE Simulation Guide},
    pdfkeywords={VLE, Molecular Dynamics, GROMACS, AmberTools, GAFF2, Methanol, Propanol},
    bookmarks=true,
    bookmarksopen=true,
    bookmarksnumbered=true
}

\begin{document}

\maketitle
\tableofcontents % Good for a long document
\newpage

\begin{abstract}
This document provides a comprehensive, step-by-step guide for simulating the Vapor-Liquid Equilibrium (VLE) of a binary organic mixture (exemplified by Methanol-Propanol) using a combination of AmberTools for parameterization and system building, and GROMACS for Molecular Dynamics (MD) simulations. It details each stage from individual molecule preparation, force field parameterization using the General Amber Force Field (GAFF2), system assembly, equilibration protocols, to VLE production runs and initial data analysis. Common pitfalls, troubleshooting advice, and explanations of key commands and concepts are included to assist users, particularly those new to MD simulations or this specific workflow. The goal is to enable the user to generate reliable VLE data points suitable for constructing phase diagrams.
\end{abstract}

\section{Introduction}
Vapor-Liquid Equilibrium (VLE) data is crucial in chemical engineering for the design and optimization of separation processes like distillation. Experimental determination of VLE data can be time-consuming and expensive. Molecular Dynamics (MD) simulations offer a powerful computational alternative to predict VLE behavior from the underlying intermolecular forces. This guide focuses on a common workflow employing the widely used GAFF2 force field (suitable for organic molecules) via AmberTools, and the efficient GROMACS MD engine.

\section{Workflow Overview}
The overall process can be visualized as follows:

\begin{figure}[h!]
\centering
\includegraphics[width=0.9\textwidth]{name (2).pdf} % Make sure 'name_2.pdf' is in your directory
% \fbox{\parbox[c][10cm][c]{0.9\textwidth}{\centering Placeholder for Workflow Image (name\_2.pdf)}}
\caption{Overall workflow for VLE simulation.}
\label{fig:workflow_overview}
\end{figure}

The flowchart depicts:
\begin{enumerate}
    \item Single Molecule Prep (PDB $\rightarrow$ antechamber $\rightarrow$ parmchk2 $\rightarrow$ .mol2 + .frcmod) - Loop for each component
    \item System Building (Packmol $\rightarrow$ .pdb\_packed $\rightarrow$ tleap $\rightarrow$ .prmtop + .inpcrd)
    \item Conversion (InterMol $\rightarrow$ .top + .gro)
    \item GROMACS Equilibration (Minimization $\rightarrow$ NVT $\rightarrow$ NPT)
    \item VLE Setup (editconf to create slab)
    \item VLE Production (NVT long run)
    \item Analysis (gmx density, gmx energy, plotting)
\end{enumerate}

\section{Software Versions}
This guide was developed and tested using the following software versions. Minor variations in newer or older versions might lead to slightly different command outputs or GUI options.

\begin{itemize}
    \item AmberTools: 22 (antechamber, parmchk2, tleap)
    \item Packmol: 20.15.1
    \item InterMol: 0.1.2 (via Conda)
    \item GROMACS: 2024.5-conda\_forge
    \item Gnuplot: 6.0 (or your preferred plotting tool)
    \item Operating System: Ubuntu Server (headless) on which simulations were run.
\end{itemize}

% --- Start of content from previous responses, correctly integrated ---
% --- Phases 1, 2, 3, 4, and 5 (up to 5.3) are assumed to be here ---

%==================================================================
%   PHASE 1: SINGLE MOLECULE PREPARATION & PARAMETERIZATION
%==================================================================
\section[Phase 1: Single Molecule Preparation \& Parameterization]{Phase 1: Single Molecule Preparation \\ \& Parameterization} % Manual line break for ToC if needed
\label{sec:phase1}
This phase must be completed for \emph{each unique chemical species} in your mixture. The goal is to obtain GAFF2-compatible parameters and a Tripos Mol2 file that describes the molecule's topology, atom types, and partial charges. For our example, this will be done for Methanol (\texttt{MOH}) and 1-Propanol (\texttt{PRH}).

\subsection{1.1. Create Initial 3D Structure (PDB format)}
\label{ssec:pdb_creation}
The very first step in parameterizing a novel molecule for molecular dynamics simulations is to obtain a reasonable three-dimensional (3D) representation of its structure. This initial structure serves as the input for \texttt{antechamber}, which will analyze its topology (which atoms are bonded to which) and geometry to assign atom types and calculate partial charges. While \texttt{antechamber} can perform some geometry adjustments, starting with a chemically sensible structure can prevent errors and lead to more reliable parameters. This structure is typically provided in the Protein Data Bank (PDB) format (\texttt{.pdb}).

\subsubsection*{Methods for Obtaining Initial PDB Structures}
There are several ways to generate a PDB file for your molecule:
\begin{itemize}
    \item \textbf{Using Molecular Sketching/Building Software (Recommended):} This is the most common and flexible approach, especially for molecules not readily available in databases or for novel structures.
    \begin{itemize}
        \item \textbf{Software Examples:}
        \begin{itemize}
            \item \textit{Avogadro:} A free, open-source, cross-platform molecular editor and visualizer. It's user-friendly and excellent for building organic molecules, adding hydrogens, and performing quick geometry clean-ups using built-in force fields like UFF or MMFF94.
            \item \textit{GaussView:} A commercial GUI for Gaussian, but also a capable molecule builder.
            \item \textit{ChemDraw/ChemDoodle:} Primarily 2D chemical drawing tools, but many versions can generate 3D coordinates and export to PDB.
            \item \textit{Jmol/BIOVIA Draw/Other Editors:} Various other free or commercial tools exist.
        \end{itemize}
        \item \textbf{General Process in a Builder (e.g., Avogadro):}
        \begin{enumerate}
            \item \textit{Build the heavy atoms:} Sketch the carbon backbone and place heteroatoms (O, N, S, etc.).
            \item \textit{Set correct bond orders:} Ensure single, double, triple, and aromatic bonds are correctly assigned. Most builders attempt this automatically but review is good.
            \item \textit{Add hydrogens:} Use the software's function to add hydrogens to satisfy valencies (e.g., Avogadro: Build $\rightarrow$ Add Hydrogens).
            \item \textit{Perform a rough geometry optimization/cleanup:} This is crucial. Use a simple, fast molecular mechanics force field (Avogadro: Extensions $\rightarrow$ Optimize Geometry). This step ensures bond lengths, angles, and dihedrals are reasonable, preventing overly strained or clashed initial structures that could cause problems for \texttt{antechamber}. This is not meant to be a high-accuracy quantum mechanical optimization.
            \item \textit{Inspect the 3D structure:} Rotate and examine the molecule to ensure it looks chemically correct.
            \item \textit{Save As PDB:} Export the structure in PDB format (\texttt{.pdb}). Pay attention to options related to atom naming or residue naming if the software provides them, though we will often standardize these manually.
        \end{enumerate}
    \end{itemize}
    \item \textbf{From Online Databases:}
    For common molecules, PDB files or similar 3D coordinate files might be available from chemical databases like:
    \begin{itemize}
        \item PubChem (\url{https://pubchem.ncbi.nlm.nih.gov/})
        \item ChemSpider (\url{http://www.chemspider.com/})
        \item NIST Chemistry WebBook (\url{https://webbook.nist.gov/chemistry/})
    \end{itemize}
    If downloading, ensure you get a 3D conformer, and verify that the structure makes sense. You might still want to open it in a molecular editor to check/assign residue and atom names as per the conventions below.
    \item \textbf{Manual Creation (Advanced Users / Very Simple Molecules):}
    It is technically possible to create a PDB file using a plain text editor. However, the PDB format has very strict column requirements for each piece of information (atom serial, atom name, residue name, coordinates, etc.). This method is highly error-prone for beginners and generally not recommended unless for extremely simple cases or by experts familiar with the format.
\end{itemize}

\subsubsection*{PDB File Content and Formatting Requirements}
A PDB file is a plain text file with specific formatting. For \texttt{antechamber} processing of a single molecule type, your PDB file should ideally contain only the atoms for one instance of that molecule.

\textbf{Key PDB Record: ATOM/HETATM lines}
These lines contain the core information. Here's a breakdown of the essential columns (1-indexed), followed by an example:
\begin{verbatim}
COLUMNS        DATA  TYPE    FIELD           DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  " or "HETATM"
 7 - 11        Integer       serial          Atom serial number.
13 - 16        Atom          name            Atom name. (e.g., " C1 ", "C1M ")
17             Character     altLoc          Alternate location indicator.
18 - 20        Residue name  resName         Residue name. (e.g., "MOH")
22             Character     chainID         Chain identifier. (e.g., "A")
23 - 26        Integer       resSeq          Residue sequence number. (e.g., "  1")
27             AChar         iCode           Code for insertion of residues.
31 - 38        Real(8.3)     x               Orthogonal coordinates for X.
39 - 46        Real(8.3)     y               Orthogonal coordinates for Y.
47 - 54        Real(8.3)     z               Orthogonal coordinates for Z.
55 - 60        Real(6.2)     occupancy       Occupancy. (e.g., "  1.00")
61 - 66        Real(6.2)     tempFactor      Temperature factor. (e.g., "  0.00")
77 - 78        LString(2)    element         Element symbol, right-justified. (e.g., " C")
79 - 80        LString(2)    charge          Charge on the atom.
\end{verbatim}
Conceptual example lines:
\begin{lstlisting}[language=tex, basicstyle=\ttfamily\footnotesize, frame=none, numbers=none]
ATOM      1  C1M MOH A   1      -1.376   0.638  -0.000  1.00  0.00           C
HETATM    5  C2P PRH A   1      -0.109  -0.201   0.000  1.00  0.00           C
\end{lstlisting}

\begin{itemize}
    \item \textbf{Cols 1-6: Record Name} (\texttt{ATOM} or \texttt{HETATM}). For general molecules, \texttt{HETATM} is often appropriate, but \texttt{ATOM} usually works too for \texttt{antechamber} input.
    \item \textbf{Cols 7-11: Atom Serial Number} (Integer, right-justified).
    \item \textbf{Cols 13-16: Atom Name} (String, typically 4 characters. Left-justified if <4 chars, e.g., \texttt{" C1 "}, \texttt{" O  "}. Spaces are significant. Our strategy of using unique 3-character names like \texttt{C1M} fits well, becoming e.g., \texttt{"C1M "}).
    \item \textbf{Cols 18-20: Residue Name} (String, 3 characters, uppercase recommended, e.g., \texttt{MOH}).
    \item \textbf{Col 22: Chain Identifier} (Character, e.g., \texttt{A}). Can be left blank for single molecules.
    \item \textbf{Cols 23-26: Residue Sequence Number} (Integer, right-justified). Usually \texttt{1} for a single, isolated molecule.
    \item \textbf{Cols 31-38: X Coordinate} (Real, F8.3 format, e.g., \texttt{-1.376}).
    \item \textbf{Cols 39-46: Y Coordinate} (Real, F8.3 format, e.g., \texttt{0.638}).
    \item \textbf{Cols 47-54: Z Coordinate} (Real, F8.3 format, e.g., \texttt{-0.000}).
    \textit{Strict adherence to these coordinate column formats is vital, as \texttt{antechamber} will fail if they are misaligned.}
    \item \textbf{Cols 55-60: Occupancy} (Real, F6.2, e.g., \texttt{1.00}). Typically \texttt{1.00}.
    \item \textbf{Cols 61-66: Temperature Factor} (Real, F6.2, e.g., \texttt{0.00}). Typically \texttt{0.00} for input structures.
    \item \textbf{Cols 77-78: Element Symbol} (String, right-justified, e.g., \texttt{" C"}, \texttt{" H"}).
\end{itemize}

\begin{showybox}{RoyalBlue}{Best Practice: Residue Naming}
    \textbf{Consistency is Key:} The residue name (e.g., \texttt{MOH} for methanol, \texttt{PRH} for propanol) assigned in the PDB file should be:
    \begin{itemize}
        \item A 3-character uppercase string (convention).
        \item Unique for each distinct molecule type in your system.
        \item The same name you will provide to \texttt{antechamber} via the \texttt{-rn} flag.
        \item The same name that will appear in the \texttt{.mol2} file generated by \texttt{antechamber}.
        \item The same name used in \texttt{tleap} when defining the residue (e.g., \texttt{MOH = loadmol2 moh.mol2}).
        \item The same name used in the \texttt{mixture\_packed.pdb} file by \texttt{packmol}.
    \end{itemize}
    \textbf{Pitfall:} Using common biochemical residue names (e.g., \texttt{ALA}, \texttt{LYS}) for non-biological molecules can sometimes cause conflicts if standard biomolecular force fields are also loaded in \texttt{tleap}. As we saw with `POL` causing issues, choosing slightly more unique names like `PRH` (Propanol Hydrocarbon-like) can prevent `tleap` from misinterpreting or "splitting" residues.
\end{showybox}

\begin{showybox}{SeaGreen}{Best Practice: Atom Naming}
    \textbf{Uniqueness within Residue:} Atom names within a single residue definition must be unique (e.g., a methanol molecule cannot have two atoms named \texttt{C1}). Standard PDB atom names are up to 4 characters.
    \textbf{Strategy for this Workflow (Global Uniqueness):}
    \begin{itemize}
        \item For this specific workflow, especially due to the Packmol PDB workaround (see Section \ref{ssec:packmol_packing}), we recommend using \textbf{globally unique atom names} across different molecule types.
        \item Example:
        \begin{itemize}
            \item Methanol (\texttt{MOH}): \texttt{C1M, H1M, H2M, H3M, O1M, HOM}
            \item Propanol (\texttt{PRH}): \texttt{C1P, H1P, H2P, C2P, H3P, H4P, C3P, H5P, H6P, H7P, O1P, HOP}
        \end{itemize}
        \item This global uniqueness (e.g., \texttt{C1M} vs \texttt{C1P}) helps \texttt{tleap} definitively assign atoms to the correct residue type when reading the \texttt{mixture\_packed.pdb} file, even if residue information were somehow ambiguous (though residue names should prevent this). It's crucial for the Packmol workaround.
    \end{itemize}
    \textbf{Formatting:} PDB atom names are typically 4 characters. If your chosen names are shorter (e.g., 3 characters like \texttt{C1M}), they should be left-justified in the 4-character field, followed by a space (e.g., \texttt{"C1M "}). Many tools handle this automatically if you specify a 3-character name.
    \textbf{Pitfall:} Inconsistent atom naming between the \texttt{.pdb} used for \texttt{antechamber}, the resulting \texttt{.mol2}, and the \texttt{.pdb} used for \texttt{packmol} will cause \texttt{tleap} to fail with "atom not found" errors.
\end{showybox}

\subsubsection*{Example: Conceptual PDB Structure for Methanol (\texttt{methanol\_unique.pdb})}
Below is a conceptual representation of what the \texttt{methanol\_unique.pdb} file might look like, emphasizing the unique atom and residue naming used in this workflow. Coordinates are illustrative.
\begin{lstlisting}[language=tex, caption=Example PDB content for Methanol (MOH), label=lst:pdb_moh, basicstyle=\ttfamily\footnotesize, numbers=none]
ATOM      1  C1M MOH A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  H1M MOH A   1       0.000   0.000   1.090  1.00  0.00           H
ATOM      3  H2M MOH A   1       1.028   0.000  -0.363  1.00  0.00           H
ATOM      4  H3M MOH A   1      -0.514  -0.890  -0.363  1.00  0.00           H
ATOM      5  O1M MOH A   1      -1.300   0.500   0.000  1.00  0.00           O
ATOM      6  HOM MOH A   1      -1.900   0.000   0.000  1.00  0.00           H
TER
END
\end{lstlisting}
This detailed PDB provides \texttt{antechamber} with a clear starting point, minimizing ambiguities in atom type assignment and connectivity.

\subsection{1.2. Generate GAFF2 Atom Types \& Charges (\texttt{antechamber})}
\label{ssec:antechamber}
Once you have a valid PDB file for your individual molecule (e.g., \texttt{methanol\_unique.pdb} with residue name \texttt{MOH}), the next step is to use the \texttt{antechamber} program from the AmberTools suite. The primary purposes of running \texttt{antechamber} at this stage are:
\begin{itemize}
    \item \textbf{Assign GAFF2 Atom Types:} \texttt{antechamber} analyzes the chemical environment of each atom (its element type, connectivity, and bond orders) and assigns it a specific GAFF2 (General Amber Force Field version 2) atom type. These atom types (e.g., \texttt{c3} for an sp3 carbon, \texttt{oh} for an alcohol oxygen, \texttt{hc} for a hydrogen bonded to an aliphatic carbon) are crucial because all subsequent force field parameters (bond lengths, angles, dihedrals, Lennard-Jones parameters) are defined based on these types.
    \item \textbf{Calculate Partial Atomic Charges:} For most simulations, especially those involving polar molecules or condensed phases, accurate representation of electrostatic interactions is vital. \texttt{antechamber} can employ various quantum mechanical (QM) methods followed by charge fitting procedures to derive a set of partial atomic charges that reproduce the molecule's electrostatic potential. A widely used and generally robust method for GAFF is AM1-BCC.
    \item \textbf{Generate a Tripos Mol2 File:} The output of this step is typically a \texttt{.mol2} file. This is a standard chemical file format that, in this context, will store the molecule's 3D coordinates, bond connectivity, the newly assigned GAFF2 atom types, the calculated partial charges, and the residue name you specified. This \texttt{.mol2} file becomes a key input for subsequent steps (\texttt{parmchk2} and \texttt{tleap}).
\end{itemize}

\subsubsection*{Command Structure}
The general command for running \texttt{antechamber} is:
\begin{lstlisting}[language=bash, caption=General antechamber command structure, numbers=none]
antechamber -i <input_file> -fi <input_format> -o <output_file> -fo <output_format> -c <charge_method> -s <verbosity> -rn <residue_name> -nc <net_charge>
\end{lstlisting}
For our Methanol (\texttt{MOH}) example:
\begin{lstlisting}[language=bash, caption=Antechamber command for Methanol (MOH), label=lst:antechamber_moh_cmd]
antechamber -i methanol_unique.pdb -fi pdb -o moh.mol2 -fo mol2 -c bcc -s 2 -rn MOH -nc 0
\end{lstlisting}

\subsubsection*{Explanation of Command Options}
\begin{description}
    \item[\texttt{-i methanol\_unique.pdb}] Specifies the input file. In our case, this is the PDB file we prepared in step 1.1.
    \item[\texttt{-fi pdb}] Specifies the format of the input file. Here, it's \texttt{pdb}.
    \item[\texttt{-o moh.mol2}] Specifies the desired output file name. It's conventional to use the molecule's identifier and the \texttt{.mol2} extension.
    \item[\texttt{-fo mol2}] Specifies the format of the output file, which is Tripos \texttt{mol2}.
    \item[\texttt{-c bcc}] Specifies the charge calculation method.
    \texttt{bcc} stands for AM1-BCC (Bond Charge Correction). This is a two-stage process: first, a semi-empirical quantum calculation (AM1) is performed to get an initial set of charges. Then, "Bond Charge Corrections" (BCCs) are applied based on bond types to refine these charges. AM1-BCC charges are generally recommended for use with the GAFF/GAFF2 force field as they tend to provide a good balance of accuracy and computational speed for a wide range of organic molecules.
    Other common options include \texttt{gas} (for Gasteiger charges, which are empirical and very fast but less accurate for condensed phases) or \texttt{resp} (Restrained Electrostatic Potential fitting, which can be more accurate but requires a QM electrostatic potential, often from a higher-level QM calculation like HF/6-31G*). For general organic molecules with GAFF2, \texttt{bcc} is a robust default.
    \item[\texttt{-s 2}] Sets the status information level (verbosity). \texttt{2} provides a good amount of information about the steps \texttt{antechamber} is taking without being overwhelming.
    \item[\texttt{-rn MOH}] Assigns the residue name in the output \texttt{.mol2} file.
    This is a critical flag. The residue name provided here (e.g., \texttt{MOH}) will be embedded in the \texttt{moh.mol2} file. Later, when \texttt{tleap} loads this \texttt{moh.mol2} file, it will define a molecular unit (residue template) with this name. It is essential for consistency, especially if your input PDB also used this residue name. We used \texttt{MOH} in our \texttt{methanol\_unique.pdb} and ensure it here too.
    \item[\texttt{-nc 0}] Specifies the net charge of the molecule.
    Methanol is a neutral molecule, so its net charge is \texttt{0}.
    For ions, you would set this to the appropriate integer charge (e.g., \texttt{1} for a sodium ion, \texttt{-1} for a chloride ion). \texttt{antechamber} will try to derive partial charges that sum up to this net charge.
\end{description}

\subsubsection*{What Happens During an \texttt{antechamber} Run?}
\begin{enumerate}
    \item \textbf{Reading Input:} Parses the input PDB file for atomic coordinates and (if present and unambiguous) connectivity.
    \item \textbf{Atom and Bond Typing (Initial):} Determines initial atom and bond types based on geometry and element.
    \item \textbf{SQM Calculation (Semi-empirical QM):} If AM1-BCC is chosen, \texttt{antechamber} calls an external semi-empirical QM program (usually \texttt{sqm}, which is part of AmberTools) to perform an AM1 calculation on the molecule. This calculation optimizes the geometry (if requested, though usually not heavily for this step) and calculates Mulliken charges and bond orders.
    You'll see files like \texttt{sqm.in}, \texttt{sqm.out}, \texttt{sqm.pdb} being created and deleted in your working directory during this step.
    \item \textbf{BCC Application:} The Mulliken charges from the AM1 calculation are then adjusted using the Bond Charge Correction parameters to yield the final AM1-BCC charges.
    \item \textbf{GAFF2 Atom Type Assignment:} Based on the (potentially QM-refined) geometry, connectivity, and element types, final GAFF2 atom types are assigned.
    \item \textbf{Writing Output:} The \texttt{moh.mol2} file is written, containing:
    \begin{itemize}
        \item The molecule name (usually taken from the residue name or input file).
        \item The number of atoms, bonds, etc.
        \item An \texttt{@<TRIPOS>ATOM} section listing each atom with its coordinates, GAFF2 atom type (e.g., \texttt{c3}, \texttt{oh}, \texttt{hc}), its residue name (e.g., \texttt{MOH}), and its calculated partial charge.
        \item An \texttt{@<TRIPOS>BOND} section listing the bonds.
    \end{itemize}
\end{enumerate}

\subsubsection*{Checking the Output (\texttt{moh.mol2})}
After \texttt{antechamber} finishes (hopefully without errors), it's good practice to open the generated \texttt{.mol2} file (e.g., \texttt{moh.mol2}) with a text editor and quickly inspect it:
\begin{itemize}
    \item Verify the molecule name at the top is what you expect (e.g., \texttt{MOH}).
    \item In the \texttt{@<TRIPOS>ATOM} section:
    \begin{itemize}
        \item Check that the assigned GAFF2 atom types (usually the 5th column after coordinates) look reasonable for your molecule. For example, an sp3 carbon should be \texttt{c3}, an alcohol oxygen \texttt{oh}, an alcohol hydrogen \texttt{ho}, hydrogens on sp3 carbons \texttt{hc} or \texttt{h1}, \texttt{h2}, etc. depending on their environment.
        \item Confirm the residue name (usually the 7th or 8th column, depends on exact Mol2 variant) is correct (e.g., \texttt{MOH}).
        \item Look at the partial charges (last column). Do they make sense chemically? (e.g., oxygen should be negative, hydroxyl hydrogen positive, carbons can be slightly positive or negative depending on neighbors). The sum of these charges should equal the \texttt{-nc} value you specified.
    \end{itemize}
\end{itemize}

\subsubsection*{Troubleshooting \texttt{antechamber}}
\begin{itemize}
    \item \textbf{PDB Formatting Errors:} As seen earlier, \texttt{antechamber} is strict about PDB coordinate columns. If it complains about this, re-check your input PDB formatting.
    \item \textbf{Unusual Valencies / Weird Geometry:} If your input PDB has chemically incorrect structures (e.g., a carbon with 5 bonds, or extremely short/long bonds), \texttt{antechamber} might fail or assign strange atom types. A quick geometry cleanup in a molecular editor often helps.
    \item \textbf{\texttt{sqm} Failures:} Sometimes the \texttt{sqm} calculation can fail for very strained or unusual molecules. The \texttt{sqm.out} file (if it remains) might contain clues. This is rarer for simple organic molecules.
    \item \textbf{Net Charge Mismatch:} If you specify an \texttt{-nc} value that is impossible for the given atoms (e.g., \texttt{-nc 1} for methane), it will likely error out or produce strange charges.
\end{itemize}
By successfully running \texttt{antechamber}, you now have a GAFF2-typed and charged representation of your molecule in a \texttt{.mol2} file, ready for the next step with \texttt{parmchk2}.

\subsection{1.3. Check for Missing Force Field Parameters (\texttt{parmchk2})}
\label{ssec:parmchk2}
After \texttt{antechamber} has successfully assigned GAFF2 atom types to your molecule and generated the \texttt{.mol2} file (e.g., \texttt{moh.mol2}), the next step is to ensure that all necessary force field parameters (for bonds, angles, and dihedrals involving these atom types) are actually present in the standard GAFF2 library.
While GAFF2 is quite comprehensive for common organic functional groups, it's possible that a specific combination of atom types in your molecule, particularly around unusual linkages or novel structures, might not have pre-defined parameters. The \texttt{parmchk2} program from AmberTools is designed to address this.

\subsubsection*{Purpose of \texttt{parmchk2}}
\begin{itemize}
    \item \textbf{Identify Missing Parameters:} \texttt{parmchk2} reads your \texttt{.mol2} file (which contains the GAFF2 atom types assigned by \texttt{antechamber}). It then compares the bonds, angles, and dihedral terms implied by your molecule's topology against the parameters available in the specified GAFF2 library.
    \item \textbf{Generate Missing Parameters by Analogy (if possible):} If parameters are missing, \texttt{parmchk2} attempts to estimate them by analogy to similar, known parameters within the GAFF2 force field. For example, if a specific C-C-O-H dihedral is missing but other similar C-C-O-X dihedrals exist, it might derive a plausible parameter.
    \item \textbf{Highlight Critically Missing Parameters:} If \texttt{parmchk2} cannot find a suitable analogy or if a parameter is deemed critical and absolutely missing, it will flag this, often with a comment like "\texttt{ATTN: NEEDS TOPOLOGY MODIFICATION}".
    \item \textbf{Create a \texttt{.frcmod} File:} The primary output of \texttt{parmchk2} is a Force Field Modification file (\texttt{.frcmod}). This file contains entries for any parameters that \texttt{parmchk2} found to be missing from the main GAFF2 library and for which it could generate or suggest values. \texttt{tleap} will later read this \texttt{.frcmod} file to supplement the main GAFF2 parameters.
\end{itemize}

\subsubsection*{Command Structure}
The general command for running \texttt{parmchk2} is:
\begin{lstlisting}[language=bash, caption=General parmchk2 command structure, numbers=none]
parmchk2 -i <input_mol2_file> -f <input_format> -o <output_frcmod_file> -s <force_field_type> [other_options]
\end{lstlisting}
For our Methanol (\texttt{MOH}) example:
\begin{lstlisting}[language=bash, caption=Parmchk2 command for Methanol (MOH), label=lst:parmchk2_moh_cmd]
parmchk2 -i moh.mol2 -f mol2 -o moh.frcmod -s gaff2
\end{lstlisting}

\subsubsection*{Explanation of Command Options}
\begin{description}
    \item[\texttt{-i moh.mol2}] Specifies the input file, which is the \texttt{.mol2} file generated by \texttt{antechamber}. This file contains the crucial GAFF2 atom types.
    \item[\texttt{-f mol2}] Specifies the format of the input file as \texttt{mol2}.
    \item[\texttt{-o moh.frcmod}] Specifies the desired output file name for the force field modifications. It's conventional to use the molecule's identifier with the \texttt{.frcmod} extension.
    \item[\texttt{-s gaff2}] Specifies the main force field type to check against. \texttt{gaff2} refers to the General Amber Force Field version 2. If you were using the original GAFF, you would use \texttt{-s gaff}.
    \item[\texttt{-a Y} or \texttt{-a N} (Optional)] Print all parameters found by \texttt{parmchk2} to standard output (\texttt{Y}) or not (\texttt{N}, default). Useful for debugging or detailed inspection but not usually needed for the \texttt{.frcmod} file itself.
    \item[\texttt{-p \$AMBERHOME/dat/leap/parm/gaff2.dat} (Optional)] Explicitly provide the path to the main GAFF2 parameter file if \texttt{parmchk2} has trouble finding it (usually not necessary if AmberTools is sourced correctly).
\end{description}

\subsubsection*{What to Expect in the Output (\texttt{moh.frcmod})}
\begin{itemize}
    \item \textbf{For common, well-parameterized molecules} (like methanol with GAFF2): The \texttt{.frcmod} file might be very short or even practically empty except for comments and section headers (\texttt{MASS, BOND, ANGLE, DIHEDRAL, IMPROPER, NONBON}). This indicates that GAFF2 already contains all the necessary parameters for the atom types found in your molecule.
    \item \textbf{For more complex or unusual molecules:} The \texttt{.frcmod} file will list any parameters that \texttt{parmchk2} generated. Each entry will specify the atom types involved and the derived parameter values (e.g., force constant, equilibrium bond length/angle, periodicity and phase for dihedrals).
    Example of a dihedral entry: \texttt{c3-c3-c3-hc 1 0.156 0.000 3.000 GAFF2, C-C-C-H, from ccnmrht (ParmBrown,1998)}
    It's good practice to review these generated parameters. If many "\texttt{ATTN:...}" messages appear, it might indicate a problem with the initial atom typing or a truly novel chemical moiety that requires more careful parameter development (beyond the scope of this automated workflow).
\end{itemize}

\subsubsection*{What \texttt{parmchk2} Does Not Do}
\begin{itemize}
    \item It does not re-calculate charges or change atom types. These come from \texttt{antechamber}.
    \item It does not guarantee the accuracy of parameters generated by analogy, especially for highly unusual systems. These are best estimates.
\end{itemize}

\subsubsection*{After Running \texttt{parmchk2}}
No direct user action is usually required with the \texttt{.frcmod} file itself, other than noting its existence. \texttt{tleap} will automatically load and use this file in a later step to ensure all parts of your molecule are described by a complete set of force field parameters. If \texttt{parmchk2} prints many "\texttt{ATTN}" messages to the screen or in the \texttt{frcmod}, it's a sign to be cautious and perhaps re-evaluate your molecule's structure or the atom types assigned by \texttt{antechamber}.
For simple molecules like methanol and propanol, \texttt{parmchk2} typically runs quickly and without issues when using GAFF2.

%==================================================================
%   END OF PHASE 1
%==================================================================

\section[Phase 2: System Building (Packing and Amber Topology)]{Phase 2: System Building \\ (Packing and Amber Topology)} % Manual line break for ToC
\label{sec:phase2}
In this phase, we will use the \texttt{.mol2} and \texttt{.frcmod} files generated for each individual component (Methanol and Propanol in our case) to construct a simulation box containing a specified number of each molecule. We will then use \texttt{tleap} to assign the complete set of GAFF2 force field parameters to this entire system.

\subsection{2.1. Pack Molecules into a Simulation Box (\texttt{packmol})}
\label{ssec:packmol_packing}

\subsubsection*{Purpose of \texttt{packmol}}
The preceding steps (Phase \ref{sec:phase1}) focused on preparing individual molecule "blueprints" – defining their structure, assigning GAFF2 atom types, calculating charges (\texttt{.mol2} files), and identifying any missing force field parameters (\texttt{.frcmod} files). However, to simulate a bulk system, especially a mixture, we need an initial three-dimensional arrangement of many copies of these molecules within a defined simulation volume (the "box").

This is where \texttt{packmol} comes in. Packmol (available at \url{http://m3g.iqm.unicamp.br/packmol}) is a powerful and widely used software package specifically designed to create initial configurations for molecular dynamics simulations. Its core function is to "pack" molecules of specified types into defined regions (like cubes, spheres, or more complex geometries) while attempting to satisfy spatial constraints, avoid severe atomic overlaps (clashes), and achieve a target number density or number of molecules.

Starting from a well-packed, reasonably dense initial configuration is crucial because:
\begin{itemize}
    \item \textbf{Reduces Equilibration Time:} A good starting point is closer to an equilibrium liquid state than, for example, a very dilute gas that would take a very long time to condense.
    \item \textbf{Prevents Extreme Forces:} If molecules are initially placed too close together, the resulting repulsive forces at the start of an energy minimization or MD run can be astronomically large, potentially causing the simulation to become unstable or "blow up." Packmol's tolerance setting helps mitigate this.
\end{itemize}

\subsubsection*{Our Specific Strategy with \texttt{packmol}}
In an ideal scenario, one might prefer to provide \texttt{packmol} with the \texttt{.mol2} files generated by \texttt{antechamber}. These \texttt{.mol2} files are rich in information, containing not only coordinates but also the GAFF2 atom types, partial charges, and bond connectivity for each defined residue (\texttt{MOH} and \texttt{PRH}). This would ensure that \texttt{packmol} is working with the most complete representation of our parameterized molecules.

However, during the development of this workflow, we encountered an issue where the specific version of \texttt{packmol} (20.15.1) available on our simulation server did not correctly process the \texttt{filetype mol2} directive in the input script. It defaulted to interpreting all input structures as PDB files, leading to errors when it tried to parse a \texttt{.mol2} file using PDB formatting rules.

To circumvent this specific \texttt{packmol} version limitation, we adopted the following robust workaround:
\begin{enumerate}
    \item \textbf{Leverage Uniquely Named PDBs:} In Phase \ref{ssec:pdb_creation}, we emphasized creating input PDB files (e.g., \texttt{methanol\_unique.pdb}, \texttt{propanol\_unique.pdb}) where each atom has a unique name that also implicitly identifies its molecule type (e.g., \texttt{C1M}, \texttt{O1M} for methanol; \texttt{C1P}, \texttt{O1P} for propanol). The residue names (\texttt{MOH}, \texttt{PRH}) were also consistently applied.
    \item \textbf{Inform \texttt{packmol} to Expect PDBs:} In the \texttt{packmol.inp} script, we explicitly set \texttt{filetype pdb}.
    \item \textbf{\texttt{packmol}'s Action:} \texttt{packmol} then reads the coordinates, atom names, and residue names directly from these \texttt{\_unique.pdb} files. It arranges the specified number of each molecule type within the defined box.
    \item \textbf{Output Consistency:} The resulting output file, \texttt{mixture\_packed.pdb}, will therefore contain all the molecules with these same unique atom names (e.g., \texttt{C1M}, \texttt{C1P}) and residue names (\texttt{MOH}, \texttt{PRH}).
    \item \textbf{\texttt{tleap} Compatibility:} This strategy is crucial because in the next phase (Section \ref{ssec:tleap_build}), \texttt{tleap} will first load our original \texttt{.mol2} files (e.g., \texttt{moh.mol2}, \texttt{pol.mol2}). These \texttt{.mol2} files were generated by \texttt{antechamber} from the \texttt{\_unique.pdb} files and thus also contain these exact same unique atom names and the correct residue definitions (\texttt{MOH}, \texttt{PRH}). When \texttt{tleap} then loads the \texttt{mixture\_packed.pdb}, it can perfectly match the atoms and residues in the packed PDB with the molecular templates it has just learned from the \texttt{.mol2} files. This ensures correct assignment of atom types, charges, and force field parameters.
\end{enumerate}
This PDB-based workaround, while a response to a specific tool behavior, highlights the general importance of maintaining consistent naming across all stages of a simulation setup pipeline. If a different, fully MOL2-compatible version of \texttt{packmol} were used, the \texttt{structure} lines in \texttt{packmol.inp} could directly point to the \texttt{.mol2} files from \texttt{antechamber}.

\subsubsection*{Creating the \texttt{packmol.inp} File}
\texttt{packmol} is controlled by a plain text input file, conventionally named with a \texttt{.inp} extension. Below is the \texttt{packmol.inp} file tailored for our methanol-propanol mixture, incorporating the PDB input strategy.

\begin{lstlisting}[language=bash, caption=Packmol input file (\texttt{packmol.inp}) for Methanol-Propanol mixture, label=lst:packmol_inp]
# Packmol input file for Methanol (MOH) and Propanol (PRH) mixture
# Version: 1.0
# Author: Mohammad Torikh
# Date: June 21, 2025
# Purpose: Generate initial packed configuration for VLE simulation.

# Global settings
# ---------------
# tolerance: Minimum distance (in Angstroms) between atoms of different molecules.
# Packmol uses an optimization algorithm to try and satisfy this.
# A value of 2.0 Angstroms is a common starting point to prevent severe initial clashes.
# Smaller values might allow for denser packing but increase clash risk.
tolerance 2.0

# filetype: Specifies the format of the coordinate files provided in the 'structure' blocks.
# Due to observed issues with 'mol2' handling in Packmol v20.15.1 on our system,
# we use 'pdb'. This necessitates that the input PDB files contain
# all necessary atom naming information consistent with what tleap will expect.
filetype pdb

# output: The name of the PDB file Packmol will generate, containing all molecules.
output mixture_packed.pdb

# Molecule definitions and packing instructions
# -------------------------------------------

# Methanol (MOH)
# 'structure' keyword indicates the start of a definition for a molecule type.
# 'methanol_unique.pdb' is the template PDB file for a single methanol molecule.
# This file must contain the residue name 'MOH' and unique atom names (e.g., C1M, O1M).
structure methanol_unique.pdb
  # 'number' specifies how many copies of this molecule to pack.
  number 200
  # 'inside cube L L L X Y Z' defines a cubic region for packing.
  # Here, '0. 0. 0.' are the coordinates of one corner of the cube.
  # '40.' is the side length of the cube in Angstroms.
  # The cube will span from (0,0,0) to (40,40,40) Angstroms.
  # The total volume is 40^3 A^3 = 64 nm^3.
  # The number of molecules and box size should be chosen to achieve a
  # density roughly appropriate for a liquid mixture. This can be estimated:
  # Mass_MOH = 200 * 32.04 g/mol; Mass_PRH = 100 * 60.10 g/mol
  # Total_mass_amu = (200*32.04 + 100*60.10) amu
  # Volume_A^3 = 40^3 A^3
  # Density ~ Total_mass_amu / Volume_A^3 * (1.66054e-24 g/amu) / (1e-24 cm^3/A^3) g/cm^3
  # For this example, it yields ~0.65 g/cm^3, a reasonable starting point for organic liquids.
  # This box will later be equilibrated under NPT conditions.
  inside cube 0. 0. 0. 40.
end structure

# Propanol (PRH)
# 'propanol_unique.pdb' is the template PDB for a single 1-propanol molecule.
# This file must contain the residue name 'PRH' and unique atom names (e.g., C1P, O1P).
structure propanol_unique.pdb
  number 100
  # Molecules of propanol will be packed into the same cubic region defined above.
  # Packmol will attempt to place all 300 molecules (200 MOH + 100 PRH)
  # within this 40x40x40 Angstrom cube.
  inside cube 0. 0. 0. 40.
end structure

# Note: Packmol offers many other options for more complex packing, such as
# different geometric constraints (spheres, boxes, around existing structures),
# fixed molecules, etc. Consult the Packmol user guide for advanced usage.
# For our VLE purpose, a simple cubic packing is sufficient as a starting point.
\end{lstlisting}

\subsubsection*{Running \texttt{packmol}}
To execute \texttt{packmol}, you redirect the input file using standard input redirection (\texttt{<}):
\begin{lstlisting}[language=bash, caption=Running Packmol, label=lst:run_packmol]
packmol < packmol.inp
\end{lstlisting}
This command should be run in the directory containing \texttt{packmol.inp}, \texttt{methanol\_unique.pdb}, and \texttt{propanol\_unique.pdb}.

\subsubsection*{Interpreting \texttt{packmol} Output}
During its run, \texttt{packmol} provides valuable feedback on the terminal:
\begin{itemize}
    \item It confirms the input file being read.
    \item It lists the types of coordinate files specified (e.g., \texttt{pdb}).
    \item It reports the output file name.
    \item It details the number of independent structures (molecule types) and their atom counts.
    \item It will show progress as it attempts to pack each molecule type and then all molecules together, often showing optimization loop information and function values related to satisfying constraints.
    \item A successful run typically ends with a "Success!" message and the final objective function value (ideally close to zero, indicating constraints were well met).
    \item It also prints a citation request for the Packmol paper.
\end{itemize}

\subsubsection*{Output Files from \texttt{packmol}}
\begin{itemize}
    \item \textbf{\texttt{mixture\_packed.pdb} (Primary Output):} This is the PDB file containing the coordinates of all packed molecules.
    \begin{itemize}
        \item It will contain \texttt{ATOM} or \texttt{HETATM} records for all atoms (in our case, 200 methanols $\times$ $\sim$6 atoms/MOH + 100 propanols $\times$ $\sim$12 atoms/PRH = $\sim$2400 atoms).
        \item Each atom record will retain the unique atom names (e.g., \texttt{C1M, O1M, C1P, O1P}) and residue names (\texttt{MOH, PRH}) from the input template PDBs.
        \item Residue sequence numbers will be assigned by \texttt{packmol} to distinguish individual molecules (e.g., the first methanol might be residue 1, the second residue 2, etc., up to 200, then propanols might start from 201).
        \item The last line of this PDB file might contain box vector information, often as a \texttt{CRYST1} record, or sometimes GROMACS-style box vectors if the Packmol version supports it well for PDB output. \texttt{tleap} will handle defining the box explicitly in the next step.
    \end{itemize}
    \item \textbf{\texttt{packmol.log} (Optional):} A log file may be generated containing more verbose details of the packing process, which can be useful for debugging if \texttt{packmol} fails or struggles to pack the system.
\end{itemize}

\subsubsection*{Quick Verification of \texttt{mixture\_packed.pdb}}
After \texttt{packmol} completes, a quick sanity check of the output PDB is recommended:
\begin{itemize}
    \item \textbf{File Size:} Ensure \texttt{mixture\_packed.pdb} has a non-zero and substantial file size.
    \item \textbf{Atom and Residue Names:}
    \begin{lstlisting}[language=bash, caption=Checking atom and residue names in packed PDB, label=lst:check_packed_pdb]
grep "ATOM " mixture_packed.pdb | grep " MOH " | head -n 5 # Check first few methanol atoms
grep "ATOM " mixture_packed.pdb | grep " PRH " | head -n 5 # Check first few propanol atoms
    \end{lstlisting}
    Confirm the residue names are \texttt{MOH} and \texttt{PRH} respectively, and that the atom names match your unique naming scheme (e.g., \texttt{C1M, O1P}).
    \item \textbf{Visual Inspection (Optional but Recommended):} If you have access to a molecular visualizer (like VMD or PyMOL) that can read PDB files, loading \texttt{mixture\_packed.pdb} can give you a visual confirmation that the molecules are packed within a box and don't have obvious, gross overlaps.
\end{itemize}
With a successful \texttt{mixture\_packed.pdb}, you have a well-defined starting configuration for your entire mixture, ready for \texttt{tleap} to build the complete molecular topology for the Amber force field.

%==================================================================
%   END OF SECTION 2.1
%==================================================================

\subsection{2.2. Generate Amber Topology \& Coordinate Files (\texttt{tleap})}
\label{ssec:tleap_build}
After successfully preparing individual molecule parameter files (\texttt{.mol2} and \texttt{.frcmod} for each component) and creating an initial packed configuration of the mixture (\texttt{mixture\_packed.pdb}), the next crucial step is to use \texttt{tleap}. \texttt{tleap} (and its graphical counterpart \texttt{xleap}) is a core program in the AmberTools suite responsible for building the complete molecular topology and assigning parameters for a given system according to the Amber force fields.

\subsubsection*{Purpose of \texttt{tleap}}
For our workflow, \texttt{tleap} will perform several key functions:
\begin{itemize}
    \item \textbf{Load Force Field:} It starts by loading a base force field, in our case, GAFF2 (\texttt{leaprc.gaff2}), which contains the standard parameters for many common atom types, bonds, angles, and dihedrals.
    \item \textbf{Load Custom Parameters:} It then loads the \texttt{.frcmod} files we generated with \texttt{parmchk2}. These files provide any missing parameters specific to our methanol (\texttt{MOH}) and propanol (\texttt{PRH}) molecules that weren't already in the base GAFF2.
    \item \textbf{Define Molecular Units (Residue Templates):} \texttt{tleap} reads our \texttt{.mol2} files (e.g., \texttt{moh.mol2}, \texttt{pol.mol2}\footnote{Note: \texttt{pol.mol2} was the filename used for propanol's \texttt{.mol2} file, which internally defined the residue as \texttt{PRH} after troubleshooting the `POL` naming issue.}). From these, it learns the definition of each of our custom molecular units (\texttt{MOH} and \texttt{PRH}), including their atom names, GAFF2 atom types, partial charges, and internal connectivity (bonds). These definitions become reusable "templates."
    \item \textbf{Assemble the System:} It reads the \texttt{mixture\_packed.pdb} file generated by \texttt{packmol}. For each atom in this PDB, \texttt{tleap} attempts to match its atom name and residue name (e.g., atom \texttt{C1M} in residue \texttt{MOH}) to the atoms within the corresponding loaded molecular unit templates (\texttt{MOH} or \texttt{PRH}).
    \item \textbf{Parameter Assignment:} Once atoms are matched, \texttt{tleap} uses the GAFF2 atom types (from the \texttt{.mol2} templates) to look up and assign all necessary force field parameters (bond stretching, angle bending, dihedral torsion, non-bonded Lennard-Jones, and electrostatic interactions using the charges from the \texttt{.mol2} files) for every interaction in the entire system.
    \item \textbf{Define Box Information:} We instruct \texttt{tleap} about the simulation box dimensions.
    \item \textbf{Output Amber System Files:} Finally, \texttt{tleap} saves two critical files:
    \begin{itemize}
        \item A parameter/topology file (\texttt{.prmtop}): This text file contains a complete description of the molecular system's topology (all atoms, residues, bonds, angles, dihedrals, exclusions) and all associated force field parameters (force constants, equilibrium values, charges, Lennard-Jones parameters). It's the "recipe book" for the simulation.
        \item An coordinate file (\texttt{.inpcrd}): This text file contains the initial atomic coordinates for all atoms in the system, usually matching those from the input PDB (\texttt{mixture\_packed.pdb}) but now in Amber's specific format and order.
    \end{itemize}
\end{itemize}

\subsubsection*{The \texttt{tleap\_mixture.in} Script}
\texttt{tleap} is controlled by a script file (conventionally with a \texttt{.in} extension). Here's the script we developed, with detailed explanations:

\begin{lstlisting}[language=leap, caption=\texttt{tleap} input script (\texttt{tleap\_mixture.in}) for Methanol-Propanol, label=lst:tleap_in, morekeywords={leaprc.gaff2, MOH, PRH, pol.mol2, moh.mol2, mixture, mixture_packed.pdb, mixture.prmtop, mixture.inpcrd, check}]
# tleap input script for creating Amber system files for Methanol-Propanol mixture
# Version: 1.0
# Author: Mohammad Torikh
# Date: June 21, 2025

# --- Load Force Field and Custom Parameters ---

# 'source' command loads another leap script or a parameter set.
# 'leaprc.gaff2' loads the General Amber Force Field version 2.
# This file defines standard atom types, bond, angle, dihedral, and LJ parameters.
# It must be sourced before loading molecule definitions that use GAFF2 atom types.
# AmberTools needs to be installed and its environment sourced for tleap to find this.
source leaprc.gaff2

# 'loadamberparams' reads a .frcmod file.
# These files contain parameters identified as missing or needing modification by parmchk2.
# We load one for each unique molecule type we parameterized.
loadamberparams moh.frcmod  # Parameters for Methanol (MOH)
loadamberparams pol.frcmod  # Parameters for Propanol (PRH).
                            # Note: This pol.frcmod was generated from pol.mol2,
                            # which itself was created by antechamber using "-rn PRH".

# --- Define Molecular Units (Residue Templates) ---

# 'loadmol2' reads a Tripos Mol2 file and uses it to define a new molecular unit (residue template).
# The variable on the left (e.g., MOH) becomes the name of this unit within tleap.
# The .mol2 file (e.g., moh.mol2) must contain:
#   1. Atom names (e.g., C1M, O1M, H1M...)
#   2. Correct GAFF2 atom types for each atom (assigned by antechamber)
#   3. Partial atomic charges for each atom (calculated by antechamber)
#   4. Bond connectivity.
#   5. The residue name embedded within the .mol2 file (from antechamber's -rn flag)
#      should ideally match the variable name here for clarity (e.g., MOH for moh.mol2).
#      When tleap loads `pol.mol2` (which was created with `antechamber ... -rn PRH`),
#      the unit defined *inside* that .mol2 file is named PRH. We assign it to a tleap
#      variable also named PRH for consistency.
MOH = loadmol2 moh.mol2
PRH = loadmol2 pol.mol2  # This pol.mol2 defines the PRH unit

# --- Load Packed Coordinates and Assemble the System ---

# 'loadpdb' reads a PDB file and creates a new system unit.
# It uses the residue and atom names from the PDB file to identify which
# pre-defined molecular units (MOH, PRH that we just loaded) to use for each
# residue encountered in the PDB.
# CRITICAL:
#   - The residue names in 'mixture_packed.pdb' (e.g., "MOH", "PRH") MUST match
#     the names of the units we defined above via loadmol2.
#   - The atom names within each residue in 'mixture_packed.pdb' (e.g., "C1M", "O1M"
#     for an MOH residue) MUST match the atom names within the corresponding
#     loaded .mol2 template (e.g., moh.mol2).
# This matching is why our unique atom naming strategy with the Packmol PDB workaround was essential.
mixture = loadpdb mixture_packed.pdb

# --- System Checks and Box Definition ---

# 'check' is a very important command. It inspects the specified unit (here, 'mixture')
# for missing parameters (bonds, angles, dihedrals) or other inconsistencies.
# Any errors or serious warnings here usually indicate a problem with atom typing,
# missing .frcmod parameters, or incorrect connectivity.
# Minor warnings about non-integral net charge for the whole system are common
# for GAFF/AM1-BCC due to floating point sums, and usually acceptable if small.
check mixture

# 'set <unit> box { X Y Z }' defines the periodic boundary box dimensions for the unit.
# The PDB file from Packmol might not contain a standard CRYST1 record that tleap uses,
# or its dimensions might not be precisely what we want to record in the prmtop.
# We set it to the dimensions used in Packmol (40.0 40.0 40.0 Angstroms).
# These dimensions will be written to the .prmtop and .inpcrd files.
# GROMACS will later read these dimensions from the converted .gro file.
set mixture box { 40.0  40.0  40.0 }

# --- Save Output Files ---

# 'saveamberparm' saves the topology (.prmtop) and coordinates (.inpcrd) for the specified unit.
# <unit>: The tleap unit to save (here, 'mixture').
# <prmtop_file>: Name of the output parameter/topology file (e.g., 'mixture.prmtop').
# <inpcrd_file>: Name of the output coordinate file (e.g., 'mixture.inpcrd').
saveamberparm mixture mixture.prmtop mixture.inpcrd

# 'quit' exits tleap.
quit
\end{lstlisting}

\subsubsection*{Running \texttt{tleap}}
\texttt{tleap} is executed by passing its input script using the \texttt{-f} flag:
\begin{lstlisting}[language=bash, caption=Running \texttt{tleap}, label=lst:run_tleap]
tleap -f tleap_mixture.in
\end{lstlisting}

\subsubsection*{Interpreting \texttt{tleap} Output and \texttt{leap.log}}
\begin{itemize}
    \item \textbf{Loading Messages:} \texttt{tleap} will report loading \texttt{leaprc.gaff2}, your \texttt{.frcmod} files, and your \texttt{.mol2} files. When loading a \texttt{.mol2} file, it will typically say "\texttt{Reading MOLECULE named PRH}" (or \texttt{MOH}), confirming the internal name of the unit defined in the \texttt{.mol2} file.
    \item \textbf{PDB Loading:} When loading \texttt{mixture\_packed.pdb}, it will indicate how many atoms it read.
    \item \textbf{\texttt{check mixture} Output:} Pay close attention to any errors or warnings here.
    \begin{itemize}
        \item "\texttt{Atom ... does not have a type.}" or "\texttt{Unknown atom type...}": This is serious. It means an atom in your system couldn't be matched to a GAFF2 type, or the type is missing parameters. This usually points to an error in your \texttt{.mol2} files or a mismatch between atom names in the PDB and the \texttt{.mol2} templates.
        \item "\texttt{Could not find bond/angle/dihedral parameter...}": This indicates missing force field parameters. Ensure your \texttt{.frcmod} files were loaded and are correct. If \texttt{parmchk2} didn't flag these as critically missing ("\texttt{ATTN}"), \texttt{tleap} might still proceed by using generic or zeroed parameters, which is undesirable.
        \item "\texttt{WARNING: The unperturbed charge of the unit (...) is not integral/zero.}": As discussed, for a large neutral system built with AM1-BCC charges, a very small, non-zero net charge (e.g., \texttt{-0.0999}) is common due to floating-point summations and usually acceptable.
    \end{itemize}
    \item \textbf{"Split Residue" Issue (Our previous problem):}
    Recall that when we initially used \texttt{POL} as the residue name for propanol, \texttt{tleap} outputted many warnings like:
    \begin{verbatim}
    Name change in pdb file residue 1 ; this residue is split into POL and OL.
    \end{verbatim}
    Followed by: "\texttt{Unknown residue: OL}" and then "\texttt{FATAL: Atom .R<OL ...>.A<HOP ...> does not have a type.}"
    This indicated that \texttt{tleap} was misinterpreting the \texttt{POL} residue name from the PDB, possibly due to a name conflict or an internal parsing rule, and was incorrectly dividing the propanol molecules. Changing the residue name to the more unique \texttt{PRH} throughout all stages (PDB input to \texttt{antechamber}, \texttt{-rn PRH} flag, \texttt{pol.mol2} content, and \texttt{mixture\_packed.pdb}) resolved this critical issue. This highlights the sensitivity of \texttt{tleap} to residue naming.
    \item \textbf{Final Messages:} A successful \texttt{tleap} run for this kind of system will typically exit with "\texttt{Errors = 0}" and a small number of warnings (often related to the net charge or unassigned chain types for non-polymeric units).
    \item \textbf{\texttt{leap.log}:} This file contains a complete transcript of the \texttt{tleap} session, including all commands, outputs, warnings, and errors. It is the primary resource for debugging \texttt{tleap} problems.
\end{itemize}

\subsubsection*{Output Files from \texttt{tleap}}
\begin{itemize}
    \item \textbf{\texttt{mixture.prmtop}:} The Amber parameter/topology file. This is a text file containing all force field parameters for the entire system.
    \item \textbf{\texttt{mixture.inpcrd}:} The Amber coordinate file, containing the initial coordinates of all atoms in the system.
    \item \textbf{\texttt{leap.log}:} The log file.
\end{itemize}
Successfully generating \texttt{mixture.prmtop} and \texttt{mixture.inpcrd} without fatal errors means you have a complete Amber-format description of your mixed system, ready for the next step: conversion to GROMACS format.

%==================================================================
%   END OF SECTION 2.2
%==================================================================

%==================================================================
%   PHASE 3: CONVERSION TO GROMACS AND SIMULATION SETUP
%==================================================================
\section[Phase 3: Conversion to GROMACS and Simulation Setup]{Phase 3: Conversion to GROMACS \\ and Simulation Setup} % Manual line break for ToC
\label{sec:phase3}
With the \texttt{mixture.prmtop} (topology and parameters) and \texttt{mixture.inpcrd} (coordinates) files successfully generated by \texttt{tleap}, we now have a complete description of our mixed molecular system in the Amber file format. However, our goal is to run the Molecular Dynamics simulations using GROMACS, which uses its own distinct file formats for topology (\texttt{.top}) and coordinates (\texttt{.gro}). Therefore, a conversion step is necessary.

\subsection{3.1. Convert Amber Files to GROMACS Format (\texttt{InterMol})}
\label{ssec:intermol_conversion}

\subsubsection*{Purpose of Conversion and Choice of InterMol}
The primary goal here is to translate the information contained in the Amber \texttt{mixture.prmtop} and \texttt{mixture.inpcrd} files into equivalent GROMACS \texttt{mixture\_intermol.top} and \texttt{mixture\_intermol.gro} files. This involves:
\begin{itemize}
    \item Mapping Amber atom types (defined in GAFF2) to GROMACS atom types.
    \item Converting bond, angle, dihedral, and improper dihedral parameters to the functional forms and units used by GROMACS.
    \item Transferring atomic charges and Lennard-Jones parameters.
    \item Writing out the coordinates in GROMACS \texttt{.gro} format.
    \item Crucially for GROMACS: Structuring the topology (\texttt{.top}) file correctly. GROMACS expects separate \texttt{[ moleculetype ]} definitions for each distinct type of molecule in the system (e.g., one for \texttt{MOH}, one for \texttt{PRH}). It then requires a \texttt{[ system ]} section to name the overall system and a \texttt{[ molecules ]} section at the end to list each \texttt{moleculetype} and the number of instances of that molecule present in the simulation.
\end{itemize}
While AmberTools' ParmEd program can perform Amber-to-GROMACS conversions using its \texttt{outparm} command, we found in our specific workflow that it tended to produce a "monolithic" GROMACS \texttt{.top} file. This means it often defined the entire system as a single large \texttt{[ moleculetype ]} rather than creating distinct \texttt{[ moleculetype ]} blocks for \texttt{MOH} and \texttt{PRH} and then listing them in a final \texttt{[ molecules ]} section. This monolithic structure leads to a "No molecules were defined in the system" error when using \texttt{gmx grompp}.

To overcome this, we turned to InterMol (\url{https://github.com/shirtsgroup/InterMol}). InterMol is a Python-based tool specifically designed for converting between various molecular simulation file formats, including Amber and GROMACS. It is generally more robust in correctly identifying and separating different residue types from an Amber \texttt{prmtop} into distinct \texttt{[ moleculetype ]} definitions in the GROMACS \texttt{.top} file and correctly generating the \texttt{[ molecules ]} section.

\subsubsection*{Installing InterMol (if not already done)}
InterMol can be installed into your Conda environment:
\begin{lstlisting}[language=bash, caption=Installing InterMol with Conda, label=lst:install_intermol]
conda install -c conda-forge intermol
\end{lstlisting}
During our setup, we confirmed that version 0.1.2 was installed via Conda.

\subsubsection*{The InterMol Python Conversion Script (\texttt{convert\_with\_intermol.py})}
InterMol can be used via its command-line interface or as a Python library. While the command-line module \texttt{python -m intermol.convert} worked reliably for our installed version (see below), for more programmatic control or if the direct CLI has issues, a Python script using its API is an alternative. The following script uses internal API components:
\begin{lstlisting}[language=python, caption=Python script for InterMol conversion (\texttt{convert\_with\_intermol.py}), label=lst:intermol_py]
# convert_with_intermol.py
# Author: Mohammad Torikh
# Date: June 21, 2025
# Purpose: Convert Amber prmtop/inpcrd to GROMACS top/gro using InterMol API.

import intermol # Main import to check version and access submodules
import intermol.amber.amber_IO as amber_IO # For reading Amber files
import intermol.gromacs.gromacs_IO as gmx_IO # For writing GROMACS files
from intermol.system import System # The central System object in InterMol

print(f"Using InterMol version: {intermol.__version__ if hasattr(intermol, '__version__') else 'Unknown (version attribute not found)'}")
print("Starting InterMol conversion...")

# Define input Amber file names
prmtop_file = 'mixture.prmtop'
coord_file = 'mixture.inpcrd'

# Define output GROMACS file names
gromacs_top_out = 'mixture_intermol.top'
gromacs_gro_out = 'mixture_intermol.gro'

# Create an empty InterMol System object. This object will be populated
# with data from the Amber files and then used to write the GROMACS files.
current_system = System(name="methanol_propanol_mixture") # Name is for internal tracking

try:
    # === Reading Amber Files ===
    print(f"Reading Amber topology from: {prmtop_file}")
    # Instantiate the Amber parameter/topology reader with the prmtop file
    amber_topology_reader = amber_IO.AmberParameters(prmtop_file)
    # Call its read method to parse the file
    amber_topology_reader.read()
    # Populate our System object with the information from the prmtop
    amber_topology_reader.fill_system(current_system)

    print(f"Reading Amber coordinates from: {coord_file}")
    # Instantiate the Amber coordinate reader with the inpcrd file
    amber_coordinate_reader = amber_IO.AmberCoordinates(coord_file)
    # Call its read method
    amber_coordinate_reader.read()
    # Populate our System object with the coordinate information
    # This links the coordinates to the topology already in current_system
    amber_coordinate_reader.fill_system(current_system)

    print("Amber files successfully loaded into InterMol System object.")

    # === Writing GROMACS Files ===
    print(f"Writing GROMACS topology to: {gromacs_top_out}")
    # Instantiate the GROMACS topology writer, passing it our populated System object
    gromacs_topology_writer = gmx_IO.GromacsTopology(current_system)
    # Call its write method, providing the desired output filename
    gromacs_topology_writer.write(gromacs_top_out)

    print(f"Writing GROMACS coordinates to: {gromacs_gro_out}")
    # Instantiate the GROMACS structure (.gro) writer
    gromacs_structure_writer = gmx_IO.GromacsStructure(current_system)
    # Call its write method
    gromacs_structure_writer.write(gromacs_gro_out)

    print("Conversion using InterMol API completed successfully.")

except FileNotFoundError as e:
    print(f"Error: Input file not found - {e}")
except ImportError as e:
    print(f"Error: Could not import an InterMol module. Is InterMol installed correctly? - {e}")
    print(f"Attempted to use 'intermol.amber.amber_IO' and 'intermol.gromacs.gromacs_IO'.")
except AttributeError as e:
    print(f"Error: An attribute was not found, this might indicate an API mismatch with your InterMol version - {e}")
except Exception as e:
    print(f"An unexpected error occurred during InterMol conversion: {e}")
\end{lstlisting}

\subsubsection*{Alternative InterMol Invocation (Command Line using \texttt{-m})}
During our troubleshooting, we found that for InterMol version 0.1.2 installed via Conda, the most straightforward way to perform the conversion was to use its module execution capability:
\begin{lstlisting}[language=bash, caption=Running InterMol via command line module, label=lst:run_intermol_cli]
python -m intermol.convert --amb_in mixture.prmtop mixture.inpcrd --gromacs --oname mixture_intermol
\end{lstlisting}
\begin{description}
    \item[\texttt{python -m intermol.convert}] This tells Python to run the \texttt{convert.py} script located within the \texttt{intermol} package.
    \item[\texttt{--amb\_in mixture.prmtop mixture.inpcrd}] Specifies the input Amber files (topology first, then coordinates).
    \item[\texttt{--gromacs}] A flag indicating that the desired output format is GROMACS.
    \item[\texttt{--oname mixture\_intermol}] Sets the base name for the output files. InterMol will create \texttt{mixture\_intermol.top} and \texttt{mixture\_intermol.gro}.
\end{description}
This command-line method is often less prone to minor API changes than writing a custom Python script that uses internal InterMol classes. For this documentation, using this command-line approach is recommended for its simplicity if it works with the installed InterMol version.

\subsubsection*{Expected Output and Verification}
InterMol should print some informational messages, including a citation request, and indicate "Finished!" upon successful completion.
The output files are:
\begin{itemize}
    \item \texttt{mixture\_intermol.top}: The GROMACS topology file.
    \item \texttt{mixture\_intermol.gro}: The GROMACS coordinate file.
\end{itemize}

\subsubsection*{Crucial Verification of \texttt{mixture\_intermol.top}}
This is the most important check. Open \texttt{mixture\_intermol.top} with a text editor or use \texttt{grep} and \texttt{tail} to verify its structure:
\begin{itemize}
    \item \textbf{Separate \texttt{[ moleculetype ]} Sections:}
    \begin{lstlisting}[language=bash, caption=Checking for moleculetype sections, numbers=none]
grep "\[ moleculetype \]" mixture_intermol.top
    \end{lstlisting}
    You \textbf{must} see this line appear twice (once for \texttt{MOH}, once for \texttt{PRH}).
    \item \textbf{Moleculetype Names:} Check the names assigned by InterMol.
    \begin{lstlisting}[language=bash, caption=Checking moleculetype names, numbers=none]
grep -A 1 "\[ moleculetype \]" mixture_intermol.top
    \end{lstlisting}
    This should show something like:
\begin{lstlisting}[language=tex, caption=Expected moleculetype names in .top file, basicstyle=\ttfamily\footnotesize, numbers=none, frame=none]
[ moleculetype ]
MOH          3 ; nrexcl
--
[ moleculetype ]
PRH          3 ; nrexcl
\end{lstlisting}
    The names \texttt{MOH} and \texttt{PRH} should match what \texttt{tleap} used.
    \item \textbf{\texttt{[ system ]} and \texttt{[ molecules ]} Sections:} Check the end of the file.
    \begin{lstlisting}[language=bash, caption=Checking system and molecules sections, numbers=none]
tail -n 10 mixture_intermol.top
    \end{lstlisting}
    You \textbf{must} find a \texttt{[ system ]} section followed by a \texttt{[ molecules ]} section that correctly lists the number of each molecule type. This section is critical for \texttt{gmx grompp}.
\begin{lstlisting}[language=tex, caption=Expected system and molecules sections in .top file, basicstyle=\ttfamily\footnotesize, numbers=none, frame=none]
[ system ]
; Name
methanol_propanol_mixture ; (or whatever name InterMol/the script assigned)

[ molecules ]
; Compound        nmols
MOH                  200
PRH                  100
\end{lstlisting}
    The names \texttt{MOH} and \texttt{PRH} here must exactly match those in the \texttt{[ moleculetype ]} definitions earlier in the file.
\end{itemize}
If these structural elements are present and correct in \texttt{mixture\_intermol.top}, then you have successfully converted your Amber system into a GROMACS-compatible format.

\subsection{3.2. Energy Minimization (GROMACS)}
\label{ssec:gmx_minimization}

\subsubsection*{Purpose and Critical Importance of Energy Minimization}
After creating an initial system configuration and converting it to the GROMACS format (\texttt{mixture\_intermol.gro} and \texttt{mixture\_intermol.top}), the system is likely to contain some unfavorable atomic arrangements:
\begin{itemize}
    \item Steric Clashes: Atoms from different molecules, or even within the same flexible molecule, might be too close together.
    \item Strained Geometries: Bond lengths, angles, or dihedrals might be distorted from their ideal equilibrium values.
    \item Unfavorable Electrostatic Interactions: Atoms with like charges might be in close proximity.
\end{itemize}
Starting a Molecular Dynamics (MD) simulation directly from such a high-energy, unrelaxed configuration would result in enormous initial forces, leading to:
\begin{itemize}
    \item Extremely large atomic displacements.
    \item Failure of constraint algorithms (like LINCS), as observed with LINCS WARNINGS and subsequent crashes in earlier attempts.
    \item Numerical instability, causing the simulation to "blow up."
\end{itemize}
Energy Minimization is a computational process that gently adjusts atomic coordinates to find a nearby local minimum on the potential energy surface. It iteratively moves atoms to reduce net forces, lowering the overall potential energy. Unlike MD, it typically doesn't involve velocities or temperature. It is an \textbf{absolutely essential step} before any equilibration or production MD run.

\subsubsection*{The \texttt{minim.mdp} File}
GROMACS simulations are controlled by parameters in an \texttt{.mdp} file. For energy minimization:
\begin{lstlisting}[language=bash, caption=GROMACS MDP file for Energy Minimization (\texttt{minim.mdp}), label=lst:minim_mdp, morekeywords={integrator, emtol, emstep, nsteps, cutoff-scheme, nstlist, rlist, coulombtype, rcoulomb, rvdw, pbc, steep, Verlet, PME, xyz}]
; minim.mdp - GROMACS MDP file for Energy Minimization
; Author: Mohammad Torikh
; Date: June 21, 2025

integrator  = steep     ; Algorithm: Steepest Descent, robust for initial minimization.
                        ; 'cg' (Conjugate Gradients) can be more efficient later.

emtol       = 1000.0    ; Tolerance for convergence (kJ/mol/nm).
                        ; Stops when max force (Fmax) < this value.
                        ; 1000.0 is common; smaller (e.g., 100.0) for more rigor.

emstep      = 0.01      ; Initial step size for minimization (nm). Max displacement per step.

nsteps      = 5000      ; Max number of minimization steps.
                        ; Stops if 'emtol' reached or 'nsteps' completed.
                        ; 5000 steps usually sufficient for Fmax < 1000 kJ/mol/nm.

; Parameters for neighbor searching, electrostatics, and VdW
; Should be generally consistent with subsequent MD runs.
cutoff-scheme   = Verlet    ; Required for modern GROMACS.
nstlist         = 20        ; Freq. to update neighbor list (every 20 steps).
rlist           = 1.2       ; Short-range neighbor list cutoff (nm).
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics.
rcoulomb        = 1.2       ; Short-range electrostatic cutoff (nm).
rvdw            = 1.2       ; Short-range van der Waals cutoff (nm).
pbc             = xyz       ; Periodic Boundary Conditions in all directions.
\end{lstlisting}

\subsubsection*{GROMACS Commands for Energy Minimization}
\begin{enumerate}
    \item \textbf{\texttt{gmx grompp} (The GROMACS Pre-processor):}
    Compiles input files (\texttt{.gro}, \texttt{.top}, \texttt{.mdp}) into a binary run input file (\texttt{.tpr}).
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{grompp} for minimization, label=lst:grompp_min]
gmx grompp -f minim.mdp -c mixture_intermol.gro -p mixture_intermol.top -o system_em.tpr -maxwarn 5
    \end{lstlisting}
    \begin{itemize}
        \item \texttt{-f minim.mdp}: Specifies the minimization parameter file.
        \item \texttt{-c mixture\_intermol.gro}: Initial GROMACS coordinate file.
        \item \texttt{-p mixture\_intermol.top}: GROMACS topology file.
        \item \texttt{-o system\_em.tpr}: Output binary run input file.
        \item \texttt{-maxwarn 5}: Allows a few non-fatal warnings.
    \end{itemize}
    \texttt{grompp} performs checks; serious issues will cause errors here.

    \item \textbf{\texttt{gmx mdrun} (The GROMACS Simulation Engine):}
    Performs the energy minimization using the \texttt{.tpr} file.
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{mdrun} for minimization, label=lst:mdrun_min]
gmx mdrun -deffnm system_em -v
    \end{lstlisting}
    \begin{itemize}
        \item \texttt{-deffnm system\_em}: Sets default base filename for input/output.
        \item \texttt{-v}: Verbose output.
    \end{itemize}
\end{enumerate}

\subsubsection*{Interpreting \texttt{mdrun} Output for Minimization}
The screen output (and \texttt{system\_em.log}) will show:
\begin{lstlisting}[language=tex, caption=Example \texttt{mdrun} output during minimization, basicstyle=\ttfamily\footnotesize, numbers=none, frame=none]
Steepest Descents:
   Tolerance (Fmax)   =  1.00000e+03
   Number of steps    =         5000
Step=    0, Dmax= 1.0e-02 nm, Epot=  7.78292e+04 Fmax= 1.79434e+05, atom= 2175
Step=    1, Dmax= 1.0e-02 nm, Epot=  7.13068e+04 Fmax= 5.88890e+04, atom= 2175
...
Step=   40, Dmax= 4.6e-03 nm, Epot=  2.49397e+03 Fmax= 9.15563e+02, atom= 2363

writing lowest energy coordinates.

Steepest Descents converged to Fmax < 1000 in 41 steps
Potential Energy  =  2.4939675e+03
Maximum force     =  9.1556299e+02 on atom 2363
Norm of force     =  8.3557272e+01
\end{lstlisting}
Key things to note:
\begin{itemize}
    \item \texttt{Epot} (Potential Energy) should generally decrease.
    \item \texttt{Fmax} (Maximum Force) should decrease significantly.
    \item Convergence message: "\texttt{Steepest Descents converged to Fmax < 1000...}" indicates success.
    \item "\texttt{writing lowest energy coordinates.}" confirms the minimized structure is saved.
\end{itemize}

\subsubsection*{Output Files from Minimization}
\begin{itemize}
    \item \textbf{\texttt{system\_em.gro} (Most Important):} Coordinates after minimization. Input for NVT.
    \item \texttt{system\_em.log}: Detailed log.
    \item \texttt{system\_em.edr}: Energy file (can plot Epot and Fmax using \texttt{gmx energy}).
    \item \texttt{system\_em.trr}: Trajectory file (if output steps were set).
\end{itemize}
If minimization fails, it might indicate severe clashes or force field issues.

\subsection{3.3. NVT Equilibration (GROMACS)}
\label{ssec:gmx_nvt}

\subsubsection*{Purpose of NVT Equilibration}
After energy minimization (\texttt{system\_em.gro}), the system is at a potential energy minimum but "cold." NVT (constant Number of particles, Volume, Temperature) equilibration aims to:
\begin{itemize}
    \item \textbf{Introduce Kinetic Energy:} Assign initial velocities corresponding to the target temperature.
    \item \textbf{Thermalize the System:} Allow kinetic energy (temperature) to distribute according to Maxwell-Boltzmann.
    \item \textbf{Allow Further Structural Relaxation:} Atoms explore conformational space at target temperature, maintaining the box volume from minimization.
\end{itemize}
This step is crucial before NPT equilibration and VLE production.

\subsubsection*{The \texttt{nvt\_equil.mdp} File}
\begin{lstlisting}[language=bash, caption=GROMACS MDP file for NVT Equilibration (\texttt{nvt\_equil.mdp}), label=lst:nvt_mdp, morekeywords={integrator, dt, nsteps, nstxout-compressed, nstvout, nstenergy, nstlog, cutoff-scheme, nstlist, rlist, coulombtype, rcoulomb, rvdw, pme_order, fourierspacing, tcoupl, tc-grps, tau_t, ref_t, pbc, gen_vel, gen_temp, gen_seed, constraints, constraint_algorithm, lincs_iter, lincs_order, md, Verlet, PME, V-rescale, xyz, h-bonds, lincs}]
; nvt_equil.mdp - GROMACS MDP file for NVT Equilibration
; Author: Mohammad Torikh
; Date: June 21, 2025
; Purpose: Bring the minimized system to the target temperature at constant volume.

integrator              = md        ; Leap-frog integrator for molecular dynamics.
dt                      = 0.002     ; Timestep in picoseconds (ps). 2 fs is standard.
nsteps                  = 50000     ; Total steps: 50000 * 0.002 ps/step = 100 ps.
                                    ; Common for initial equilibration: 100 ps to 1 ns.

; Output control (nst* values are in number of steps)
nstxout-compressed      = 500       ; Save compressed coordinates (.xtc) every 1 ps.
nstvout                 = 0         ; Do not save velocities for equilibration.
nstenergy               = 500       ; Save energy data to .edr every 1 ps.
nstlog                  = 500       ; Write to log file (.log) every 1 ps.

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20        ; Freq. to update neighbor list (20 steps = 40 fs).
rlist                   = 1.2       ; Short-range neighbor list cutoff (nm).

; Electrostatics and Van der Waals
coulombtype             = PME
rcoulomb                = 1.2       ; Short-range electrostatic cutoff (nm).
rvdw                    = 1.2       ; Short-range Van der Waals cutoff (nm).
pme_order               = 4         ; Interpolation order for PME (cubic).
fourierspacing          = 0.16      ; Grid spacing for PME FFT (nm).

; Temperature coupling (Thermostat) - CRUCIAL for NVT
tcoupl                  = V-rescale ; Velocity rescaling thermostat. Robust.
tc-grps                 = System    ; Couple all atoms.
tau_t                   = 0.1       ; Time constant for temperature coupling (ps).
ref_t                   = 298.15    ; Reference temperature (Kelvin).

; Pressure coupling (Barostat) - NOT USED IN NVT
; pcoupl                = No        ; Explicitly turn off.

; Periodic Boundary Conditions
pbc                     = xyz

; Initial Velocity Generation
gen_vel                 = yes       ; Generate velocities from Maxwell distribution.
gen_temp                = 298.15    ; Temperature for Maxwell distribution (K), match ref_t.
gen_seed                = -1        ; Seed for random number generator (-1 for unique).

; Constraints - Essential for 2 fs timestep
constraints             = h-bonds   ; Constrain bonds involving hydrogen.
constraint_algorithm    = lincs     ; LINCS algorithm.
lincs_iter              = 1         ; Number of iterations for LINCS.
lincs_order             = 4         ; Order of LINCS.
\end{lstlisting}

\subsubsection*{GROMACS Commands for NVT Equilibration}
\begin{enumerate}
    \item \textbf{\texttt{gmx grompp}:}
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{grompp} for NVT equilibration, label=lst:grompp_nvt]
gmx grompp -f nvt_equil.mdp -c system_em.gro -p mixture_intermol.top -o system_nvt.tpr -maxwarn 5
    \end{lstlisting}
    Uses \texttt{system\_em.gro} (minimized coordinates).

    \item \textbf{\texttt{gmx mdrun}:}
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{mdrun} for NVT equilibration, label=lst:mdrun_nvt]
gmx mdrun -deffnm system_nvt -v
    \end{lstlisting}
\end{enumerate}

\subsubsection*{Monitoring NVT Equilibration}
\begin{itemize}
    \item \textbf{Temperature:} Use \texttt{gmx energy -f system\_nvt.edr -o temperature\_nvt.xvg} (select Temperature). Should stabilize around \texttt{ref\_t}.
    \item \textbf{Potential Energy:} Use \texttt{gmx energy -f system\_nvt.edr -o potential\_energy\_nvt.xvg} (select Potential). Should reach a stable fluctuation.
    \item \textbf{LINCS Warnings:} Check \texttt{system\_nvt.log}. Should be absent or rare after minimization.
\end{itemize}

\subsubsection*{Output Files from NVT Equilibration}
\begin{itemize}
    \item \texttt{system\_nvt.gro}: Coordinates at the end of NVT.
    \item \texttt{system\_nvt.xtc}: Trajectory.
    \item \texttt{system\_nvt.edr}: Energy file.
    \item \textbf{\texttt{system\_nvt.cpt} (Most Important for next step):} Checkpoint file, saves full system state.
\end{itemize}
Once temperature and potential energy are stable, the system is ready for NPT.

\subsection{3.4. NPT Equilibration (GROMACS)}
\label{ssec:gmx_npt}

\subsubsection*{Purpose of NPT Equilibration}
Following NVT, NPT (constant Number of particles, Pressure, Temperature) equilibration aims to:
\begin{itemize}
    \item \textbf{Achieve Target Pressure and Density:} Allow system volume to fluctuate to reach equilibrium density at the specified target pressure and temperature.
    \item \textbf{Further System Relaxation:} Box size changes allow more extensive structural relaxation.
    \item \textbf{Prepare for Slab Creation:} Critical to start VLE with a liquid phase at the correct, stable density.
\end{itemize}

\subsubsection*{The \texttt{npt.mdp} File}
Key consideration: choice of barostat. Berendsen is good for initial stability; Parrinello-Rahman for production NPT if needed later. We use Berendsen here.
\begin{lstlisting}[language=bash, caption=GROMACS MDP file for NPT Equilibration (\texttt{npt.mdp}), label=lst:npt_mdp, morekeywords={integrator, dt, nsteps, nstxout-compressed, nstenergy, nstlog, cutoff-scheme, nstlist, rlist, coulombtype, rcoulomb, rvdw, pme_order, fourierspacing, tcoupl, tc-grps, tau_t, ref_t, pcoupl, pcoupltype, tau_p, ref_p, compressibility, pbc, gen_vel, constraints, constraint_algorithm, lincs_iter, lincs_order, md, Verlet, PME, V-rescale, Berendsen, isotropic, xyz, h-bonds, lincs}]
; npt.mdp - GROMACS MDP file for NPT Equilibration
; Author: Mohammad Torikh
; Date: June 21, 2025
; Purpose: Equilibrate system density at target temperature and pressure.

integrator              = md
dt                      = 0.002     ; 2 fs timestep
nsteps                  = 500000    ; 1 ns (500,000 * 0.002 ps). Monitor density; may need 1-10 ns.

; Output control
nstxout-compressed      = 5000      ; Save coordinates every 10 ps
nstenergy               = 5000      ; Save energies every 10 ps
nstlog                  = 5000      ; Update log file every 10 ps

; Neighbor searching (can be same as NVT)
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2

; Electrostatics and VdW (can be same as NVT)
coulombtype             = PME
rcoulomb                = 1.2
rvdw                    = 1.2
pme_order               = 4
fourierspacing          = 0.16

; Temperature coupling (Thermostat)
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 298.15    ; Target temperature (Kelvin)

; Pressure coupling (Barostat) - CRUCIAL for NPT
pcoupl                  = Berendsen ; Berendsen barostat for stable equilibration.
pcoupltype              = isotropic ; Couple pressure isotropically.
tau_p                   = 2.0       ; Time constant for pressure coupling (ps). Typical: 0.5-5 ps.
ref_p                   = 1.0       ; Reference pressure (bar).
compressibility         = 4.5e-5    ; Compressibility of system (bar^-1). Water's value is often used.

; Periodic Boundary Conditions
pbc                     = xyz

; Initial Velocity Generation
gen_vel                 = no        ; IMPORTANT: Continue from NVT checkpoint file (.cpt).

; Constraints
constraints             = h-bonds
constraint_algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4
\end{lstlisting}

\subsubsection*{GROMACS Commands for NPT Equilibration}
\begin{enumerate}
    \item \textbf{\texttt{gmx grompp}:}
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{grompp} for NPT equilibration, label=lst:grompp_npt]
gmx grompp -f npt.mdp -c system_nvt.gro -p mixture_intermol.top -t system_nvt.cpt -o system_npt.tpr -maxwarn 5
    \end{lstlisting}
    \begin{itemize}
        \item \texttt{-c system\_nvt.gro}: Starting coordinates from NVT.
        \item \texttt{-t system\_nvt.cpt}: Crucial for continuity; reads full state from NVT.
    \end{itemize}
    \item \textbf{\texttt{gmx mdrun}:}
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{mdrun} for NPT equilibration, label=lst:mdrun_npt]
gmx mdrun -deffnm system_npt -v
    \end{lstlisting}
\end{enumerate}

\subsubsection*{Monitoring NPT Equilibration}
\begin{itemize}
    \item \textbf{Density:} Most important. Use \texttt{gmx energy -f system\_npt.edr -o density\_npt.xvg} (select Density). Must converge to a stable average.
    \item \textbf{Pressure:} Check with \texttt{gmx energy -f system\_npt.edr -o pressure\_npt.xvg} (select Pressure). Should fluctuate around \texttt{ref\_p}.
    \item \textbf{Potential Energy, Temperature, Volume:} Also monitor for stability.
\end{itemize}
\textbf{Extending NPT Runs (if needed):} If density isn't stable, extend the run.
\begin{enumerate}
    \item Update \texttt{npt.mdp} for more steps.
    \item Use \texttt{gmx convert-tpr} to create a new \texttt{.tpr} from the previous one:
    \begin{lstlisting}[language=bash, caption=Extending NPT run with \texttt{convert-tpr}, label=lst:extend_npt_converttpr]
# Example: extend by another 1 ns (500,000 steps if dt=0.002)
gmx convert-tpr -s system_npt.tpr -extend 1000 -o system_npt_extend.tpr
# -extend value is in ps
    \end{lstlisting}
    \item Run \texttt{mdrun} with the new \texttt{.tpr}, continuing from the previous checkpoint:
    \begin{lstlisting}[language=bash, caption=Running extended NPT simulation, label=lst:run_extended_npt]
gmx mdrun -deffnm system_npt_extend -cpi system_npt.cpt -v
    \end{lstlisting}
    Alternatively, use \texttt{gmx mdrun -deffnm system\_npt -cpi system\_npt.cpt -append} if \texttt{nsteps} in the original \texttt{.tpr} was large enough or \texttt{mdrun} is allowed to run longer by other means.
\end{enumerate}

\subsubsection*{Output Files from NPT Equilibration}
\begin{itemize}
    \item \texttt{system\_npt.gro}: Coordinates at end of NPT (should have equilibrated density).
    \item \texttt{system\_npt.xtc}: Trajectory.
    \item \texttt{system\_npt.edr}: Energy file (density, pressure, volume, etc.).
    \item \texttt{system\_npt.cpt}: Final checkpoint file.
\end{itemize}
Once density is stable, the system is ready for slab creation.

% --- END OF PHASE 3 --- (Assumes previous content for Phases 1, 2, and 3 is in place)

%==================================================================
%   PHASE 4: VLE SIMULATION SETUP & PRODUCTION
%==================================================================
\section[Phase 4: VLE Simulation Setup \& Production]{Phase 4: VLE Simulation Setup \\ \& Production} % Manual line break for ToC
\label{sec:phase4}
Having successfully equilibrated our binary mixture in the NPT ensemble (resulting in \texttt{system\_npt.gro} with an equilibrated liquid density), we are now ready to create the specific geometry required to observe VLE and then run the long production simulation.

\subsection{4.1. Create Slab Geometry (\texttt{gmx editconf})}
\label{ssec:gmx_slab}

\subsubsection*{Purpose of the Slab Geometry for VLE Simulations}
As discussed previously, to simulate the coexistence of liquid and vapor phases, we need to create a system where both can physically exist and interact. The "slab geometry" is the standard approach for this in MD simulations. It involves:
\begin{itemize}
    \item Starting with a well-equilibrated bulk liquid phase (from our NPT run).
    \item Placing this liquid slab into a simulation box that is significantly elongated in one dimension (typically chosen as the Z-axis).
    \item This elongation creates empty regions (vacuum) above and below the liquid slab along the Z-axis.
\end{itemize}
When an MD simulation is run on this slab system (usually under NVT conditions, as the overall box volume is now fixed and large), molecules from the liquid surface will evaporate into the vacuum regions, eventually forming a distinct vapor phase. Over sufficient simulation time, a dynamic equilibrium will be established, with molecules continuously transitioning between the liquid and vapor phases. This setup allows us to directly sample the properties of both coexisting phases.

\subsubsection*{The \texttt{gmx editconf} Tool}
GROMACS provides the \texttt{gmx editconf} tool to manipulate simulation box dimensions and molecule orientations. We will use it to take our NPT-equilibrated box and stretch it along the Z-axis.

\subsubsection*{Steps to Create the Slab}
\begin{enumerate}
    \item \textbf{Determine the Final Box Dimensions from the NPT Run:}
    The last line of your \texttt{system\_npt.gro} file contains the final box vectors from the NPT equilibration. For an orthogonal box, these are the X, Y, and Z dimensions of the liquid phase at equilibrium.
    \begin{lstlisting}[language=bash, caption=Getting final box dimensions from \texttt{system\_npt.gro}, label=lst:tail_npt_box]
tail -n 1 system_npt.gro
    \end{lstlisting}
    Let's assume this command outputs, for example: \texttt{3.03201  3.03201  3.03201}
    So, $X_{\text{liquid}} = 3.03201$~nm, $Y_{\text{liquid}} = 3.03201$~nm, and $Z_{\text{liquid}} = 3.03201$~nm.

    \item \textbf{Decide on the New Elongated Z Dimension ($Z_{\text{new}}$):}
    The Z-dimension of the slab box ($Z_{\text{new}}$) needs to be large enough to accommodate the liquid slab and two sufficiently large vapor regions to minimize finite-size effects or interactions between the periodic images of the slab across the Z-boundary.
    A common rule of thumb is to make $Z_{\text{new}}$ at least 3 to 5 times $Z_{\text{liquid}}$.
    For our example, if $Z_{\text{liquid}} = 3.03201$~nm, let's choose $Z_{\text{new}} = 4 \times Z_{\text{liquid}}$.
    $Z_{\text{new}} = 4 \times 3.03201 \text{ nm} = 12.12804 \text{ nm}$.
    The X and Y dimensions will remain the same as $X_{\text{liquid}}$ and $Y_{\text{liquid}}$.

    \item \textbf{Calculate Centering Coordinates for the New Box:}
    We want our liquid slab to be centered within this new, elongated box. The centering coordinates are half of the new box dimensions.
    center\_X = $X_{\text{liquid}} / 2 = 3.03201 / 2 = 1.516005$~nm
    center\_Y = $Y_{\text{liquid}} / 2 = 3.03201 / 2 = 1.516005$~nm
    center\_Z$_{\text{new}}$ = $Z_{\text{new}} / 2 = 12.12804 / 2 = 6.06402$~nm

    \item \textbf{Execute \texttt{gmx editconf}:}
    \begin{lstlisting}[language=bash, caption=Using \texttt{gmx editconf} to create slab geometry, label=lst:gmx_editconf_slab]
gmx editconf -f system_npt.gro -o system_slab_initial.gro -box 3.03201 3.03201 12.12804 -center 1.516005 1.516005 6.06402
    \end{lstlisting}
    \textbf{Explanation of Options:}
    \begin{itemize}
        \item \texttt{-f system\_npt.gro}: Specifies the input coordinate file (final structure from NPT).
        \item \texttt{-o system\_slab\_initial.gro}: Output coordinate file with slab geometry.
        \item \texttt{-box $X_{\text{liquid}}$ $Y_{\text{liquid}}$ $Z_{\text{new}}$}: Sets new box dimensions (nm).
        \item \texttt{-center center\_X center\_Y center\_Z$_{\text{new}}$}: Places the geometric center of molecules at these coordinates in the new box, ensuring the slab is centered.
    \end{itemize}
\end{enumerate}

\subsubsection*{Output of \texttt{gmx editconf}}
The command will print information about the old box, the shift applied, and the new box dimensions/volume. Verify that the "new box vectors" match your specifications.
Example output:
\begin{lstlisting}[language=tex, caption=Example output from \texttt{gmx editconf}, basicstyle=\ttfamily\footnotesize, numbers=none, frame=none]
new center      :  1.516  1.516  6.064 (nm)
new box vectors :  3.032  3.032 12.128 (nm)
new box volume  : 111.49               (nm^3)
\end{lstlisting}
The primary output is \texttt{system\_slab\_initial.gro}. This \texttt{.gro} file now contains your molecules arranged as a liquid slab, centered in a box elongated in Z.

\subsubsection*{Verification (Recommended)}
\begin{itemize}
    \item \textbf{Check Box Dimensions:}
    \begin{lstlisting}[language=bash, caption=Checking new box dimensions from \texttt{system\_slab\_initial.gro}, label=lst:tail_slab_box]
tail -n 1 system_slab_initial.gro
    \end{lstlisting}
    The last line should reflect the new elongated box dimensions (e.g., \texttt{3.03201  3.03201 12.12804}).
    \item \textbf{Visual Inspection:} Load \texttt{system\_slab\_initial.gro} into a molecular visualizer (e.g., VMD, PyMOL). You should clearly see:
    \begin{itemize}
        \item Compact X and Y dimensions.
        \item A much longer Z dimension.
        \item Molecules clustered as a "slab" in the middle of the Z-range.
        \item Significant empty (vacuum) space above and below the slab along Z.
    \end{itemize}
\end{itemize}
This \texttt{system\_slab\_initial.gro} file is the starting point for your VLE production simulation.

\subsection{4.2. VLE Production Run (NVT, GROMACS)}
\label{ssec:gmx_vle_prod}

\subsubsection*{Purpose of the VLE Production Run}
With the slab geometry established (\texttt{system\_slab\_initial.gro}), this long Molecular Dynamics simulation is where the actual Vapor-Liquid Equilibrium is simulated and data for its characterization is collected. Key objectives:
\begin{itemize}
    \item \textbf{Phase Separation:} Molecules from the liquid slab evaporate into vacuum regions, forming a vapor phase, while vapor molecules condense back.
    \item \textbf{Equilibrium Establishment:} Over time, evaporation and condensation rates equalize, leading to dynamic equilibrium. Properties of coexisting phases (densities, compositions, pressure) stabilize.
    \item \textbf{Data Collection:} Generate trajectory (\texttt{.xtc} or \texttt{.trr}) and energy (\texttt{.edr}) files for analysis in Phase \ref{sec:phase5}.
\end{itemize}
This simulation is typically run under NVT conditions (constant Number of particles, Volume, Temperature).
\begin{itemize}
    \item \textbf{N (Number of particles):} Constant.
    \item \textbf{V (Volume):} Overall simulation box volume (from \texttt{gmx editconf}) is constant.
    \item \textbf{T (Temperature):} A thermostat maintains the desired VLE temperature.
\end{itemize}
Pressure is not explicitly controlled. The pressure that develops, particularly in the equilibrated vapor phase, is the vapor pressure ($P_{\text{vap}}$) of the mixture at the chosen temperature.

\subsubsection*{The \texttt{vle\_nvt.mdp} File}
Similar to \texttt{nvt\_equil.mdp} but with a significantly longer run time (\texttt{nsteps}) and potentially adjusted output frequencies.
\begin{lstlisting}[language=bash, caption=GROMACS MDP file for VLE Production Run (\texttt{vle\_nvt.mdp}), label=lst:vle_nvt_mdp, morekeywords={integrator, dt, nsteps, nstxout-compressed, nstvout, nstfout, nstenergy, nstlog, cutoff-scheme, nstlist, rlist, coulombtype, rcoulomb, rvdw, pme_order, fourierspacing, fourier_nx, fourier_ny, fourier_nz, tcoupl, tc-grps, tau_t, ref_t, pbc, gen_vel, gen_temp, gen_seed, constraints, constraint_algorithm, lincs_iter, lincs_order, md, Verlet, PME, V-rescale, xyz, h-bonds, lincs}]
; vle_nvt.mdp - GROMACS MDP file for VLE Production Run
; Author: Mohammad Torikh
; Date: June 21, 2025
; Purpose: Simulate two-phase coexistence to determine VLE properties.

integrator              = md        ; Leap-frog integrator.
dt                      = 0.002     ; Timestep (ps) (2 fs).
nsteps                  = 10000000  ; Total steps: 10,000,000 * 0.002 ps/step = 20 ns.
                                    ; VLE often requires 10s to 100s of ns.
                                    ; For testing, use shorter run (e.g., 1-2 ns).

; Output control - Write less frequently for long runs.
nstxout-compressed      = 10000     ; Save .xtc every 20 ps (10000 steps).
nstvout                 = 0         ; No velocities for full VLE trajectory.
nstfout                 = 0         ; No forces.
nstenergy               = 10000     ; Save .edr every 20 ps.
nstlog                  = 10000     ; Write to .log every 20 ps.

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20        ; Update neighbor list every 40 fs.
rlist                   = 1.2       ; Short-range neighbor list cutoff (nm).

; Electrostatics and Van der Waals
coulombtype             = PME
rcoulomb                = 1.2
rvdw                    = 1.2
pme_order               = 4
fourierspacing          = 0.16      ; GROMACS determines PME grid. For anisotropic boxes
                                    ; (like slab), grompp might choose e.g. 20x20x80.
                                    ; If PME load is high (>35-40%), manual tuning of
                                    ; fourier_nx, fourier_ny, fourier_nz or fourierspacing
                                    ; might be needed for optimization.

; Temperature coupling (Thermostat)
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1       ; Coupling time constant (ps).
ref_t                   = 298.15    ; Target temperature (K) for VLE. Key variable.

; Pressure coupling - NONE for NVT
; pcoupl                  = No

; Periodic Boundary Conditions
pbc                     = xyz       ; 3D PBC essential.

; Initial Velocity Generation
gen_vel                 = yes       ; Generate velocities from Maxwell distribution.
                                    ; Starting from system_slab_initial.gro (no velocities).
gen_temp                = 298.15    ; Temperature for Maxwell distribution (K), match ref_t.
gen_seed                = -1        ; Random seed. Use integer for reproducibility.

; Constraints
constraints             = h-bonds
constraint_algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4
\end{lstlisting}

\subsubsection*{Explanation of Key \texttt{vle\_nvt.mdp} Parameters}
\begin{itemize}
    \item \textbf{\texttt{nsteps} (Critical):} Determines total simulation time. VLE often requires long timescales (tens to hundreds of ns). A 20 ns run (\texttt{nsteps = 10000000} with \texttt{dt = 0.002}) is a reasonable start, but longer might be needed. For initial testing, use a much shorter \texttt{nsteps}.
    \item \textbf{\texttt{nstxout-compressed, nstenergy, nstlog}:} Output frequencies reduced for long runs to save disk space (e.g., every 10-50 ps).
    \item \textbf{\texttt{fourierspacing} / PME Grid:} \texttt{grompp} sets PME grid. For anisotropic boxes, PME can be a bottleneck. If "Estimate for the relative computational load of the PME mesh part" from \texttt{grompp} is high ($>$0.4), manual tuning of \texttt{fourier\_nx, fourier\_ny, fourier\_nz} or \texttt{fourierspacing} might be needed for better load balancing in future optimized runs.
    \item \textbf{\texttt{ref\_t}:} Temperature for VLE properties. To build a T-xy diagram, repeat VLE simulation at different \texttt{ref\_t} values.
    \item \textbf{\texttt{gen\_vel = yes}:} Since \texttt{system\_slab\_initial.gro} is static, GROMACS generates initial velocities.
\end{itemize}

\subsubsection*{GROMACS Commands for VLE Production Run}
\begin{enumerate}
    \item \textbf{\texttt{gmx grompp}:}
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{grompp} for VLE production, label=lst:grompp_vle]
gmx grompp -f vle_nvt.mdp -c system_slab_initial.gro -p mixture_intermol.top -o system_vle.tpr -maxwarn 5
    \end{lstlisting}
    Uses \texttt{system\_slab\_initial.gro} (slab coordinates) and \texttt{mixture\_intermol.top}.
    Review \texttt{grompp} output, especially PME grid/load estimate.

    \item \textbf{\texttt{gmx mdrun} (The Long Run):}
    \begin{lstlisting}[language=bash, caption=GROMACS \texttt{mdrun} for VLE production, label=lst:mdrun_vle]
nohup gmx mdrun -deffnm system_vle -v > mdrun_vle.out 2>&1 &
    \end{lstlisting}
    \begin{itemize}
        \item \texttt{-deffnm system\_vle}: Uses \texttt{system\_vle.tpr}; produces \texttt{system\_vle.log, .edr, .xtc, .gro, .cpt}.
        \item \texttt{nohup ... \&}: Recommended for background execution on a server, redirecting output.
    \end{itemize}
\end{enumerate}

\subsubsection*{Monitoring the VLE Production Run}
\begin{itemize}
    \item \textbf{\texttt{mdrun\_vle.out} and \texttt{system\_vle.log}:} Use \texttt{tail -f} to monitor progress.
    \item \textbf{Trajectory Visualization (Occasionally):} After significant time (e.g., few ns), transfer \texttt{system\_vle.xtc} and a coordinate/topology file (\texttt{system\_slab\_initial.gro} or \texttt{system\_vle.tpr}) to visualize with VMD/PyMOL. Check for:
    \begin{itemize}
        \item Clear interface formation.
        \item Molecules in the vapor phase.
        \item System stability (no blow-ups, liquid not drifting due to PBC issues if X/Y box too small).
    \end{itemize}
    \item \textbf{Key Energies/Properties (using \texttt{gmx energy}):}
    \begin{itemize}
        \item Temperature: Should fluctuate around \texttt{ref\_t}.
        \item Pressure: Will fluctuate. Average over equilibrated vapor phase is $P_{\text{vap}}$.
        \item Potential Energy and Density: Will change during interface formation, then stabilize. Density will show two plateaus in the profile.
    \end{itemize}
\end{itemize}
This VLE production run is where the system evolves towards equilibrium. Patience is key for these computationally intensive runs.

% --- END OF PHASE 4 --- (Assumes previous content for Phases 1, 2, 3, and 4 is in place)

%==================================================================
%   PHASE 5: ANALYSIS AND INTERPRETATION
%==================================================================
\section[Phase 5: Analysis and Interpretation]{Phase 5: Analysis and Interpretation} % Manual line break for ToC
\label{sec:phase5}
After your long VLE production run (\texttt{gmx mdrun -deffnm system\_vle}) has completed successfully, you will have a set of output files (primarily \texttt{system\_vle.xtc}, \texttt{system\_vle.edr}, and \texttt{system\_vle.tpr}) that contain the raw data about your system's behavior at the simulated temperature. The goal of this phase is to analyze this data to determine the properties of the coexisting liquid and vapor phases, which constitute one point on your VLE diagram.

\subsection{5.1. Density Profiles (\texttt{gmx density})}
\label{ssec:gmx_density_profiles}

\subsubsection*{Purpose}
The most direct way to confirm the formation of distinct liquid and vapor phases and to determine their boundaries and compositions is by calculating the density profile of each component (and the total system) along the Z-axis (the dimension perpendicular to the liquid-vapor interfaces).
\begin{itemize}
    \item A high-density plateau in the center of the box indicates the liquid slab.
    \item Low-density plateaus at the ends of the box indicate the vapor phases.
    \item The regions of transition between these plateaus are the liquid-vapor interfaces.
\end{itemize}

\subsubsection*{The \texttt{gmx density} Tool}
GROMACS provides the \texttt{gmx density} tool for this purpose.

\subsubsection*{Commands to Generate Density Profiles}
You will run \texttt{gmx density} separately for each component and for the total system. It's important to analyze a significant portion of your trajectory, potentially discarding an initial part if the system was still evolving towards phase equilibrium (e.g., if interface formation took some nanoseconds). The \texttt{-b} (begin time) and \texttt{-e} (end time) flags in \texttt{gmx density} can be used for this, but for a first pass, analyzing the whole trajectory is common.

\textbf{For Methanol (\texttt{MOH}):}
\begin{lstlisting}[language=bash, caption=Generating methanol density profile with \texttt{gmx density}, label=lst:gmx_density_moh]
gmx density -f system_vle.xtc -s system_vle.tpr -d Z -o density_Z_MOH.xvg -sl 120
\end{lstlisting}
\begin{itemize}
    \item \texttt{-f system\_vle.xtc}: Input trajectory file.
    \item \texttt{-s system\_vle.tpr}: Input run file (topology and system information).
    \item \texttt{-d Z}: Calculate density along the Z-axis.
    \item \texttt{-o density\_Z\_MOH.xvg}: Output file for methanol density profile (\texttt{.xvg} is Grace/XMGRACE format).
    \item \texttt{-sl 120}: Number of slices along Z-axis. For a $\approx 12$~nm box, 120 slices give 0.1~nm resolution. Adjust for reasonable resolution (0.05 to 0.1~nm/slice often good).
\end{itemize}
\textbf{Interactive Prompt:} \texttt{gmx density} will list atom groups (e.g., System, Other, MOH, PRH). Select the group for Methanol (e.g., by typing its number). In our case, \texttt{MOH} was group 2.

\textbf{For Propanol (\texttt{PRH}):}
\begin{lstlisting}[language=bash, caption=Generating propanol density profile with \texttt{gmx density}, label=lst:gmx_density_prh]
gmx density -f system_vle.xtc -s system_vle.tpr -d Z -o density_Z_PRH.xvg -sl 120
\end{lstlisting}
When prompted, select the group for Propanol (e.g., \texttt{PRH}, group 3).

\textbf{For the Total System:}
\begin{lstlisting}[language=bash, caption=Generating total system density profile with \texttt{gmx density}, label=lst:gmx_density_total]
gmx density -f system_vle.xtc -s system_vle.tpr -d Z -o density_Z_total.xvg -sl 120
\end{lstlisting}
When prompted, select the group for the entire System (usually group 0 or 1).

\subsubsection*{Output \texttt{.xvg} Files}
Each command produces an \texttt{.xvg} file (plain text):
\begin{itemize}
    \item Lines starting with \texttt{\#} or \texttt{@} are comments/metadata.
    \item Data lines:
    \begin{itemize}
        \item Column 1: Z-coordinate (nm) of the slice center.
        \item Column 2: Density (kg/m$^3$) in that slice, averaged over analyzed frames.
    \end{itemize}
\end{itemize}

\subsubsection*{Visualizing Density Profiles (e.g., with Gnuplot)}
Use Gnuplot to visualize profiles and save images.
\begin{lstlisting}[language=gnuplot, caption=Gnuplot script for visualizing density profiles, label=lst:gnuplot_density]
# Start gnuplot:
# gnuplot

# Plotting individual profiles and saving:
set terminal pngcairo size 800,600 enhanced font "arial,10"
set output 'density_profile_MOH.png'
set title "Density Profile of Methanol (MOH) along Z-axis"
set xlabel "Z-coordinate (nm)"
set ylabel "Density (kg/m^3)"
set grid
set yrange [0:*] # Ensure y-axis starts at 0, auto-scale max
plot 'density_Z_MOH.xvg' using 1:2 with lines linewidth 2 title 'Methanol'
unset output

set output 'density_profile_PRH.png'
set title "Density Profile of Propanol (PRH) along Z-axis"
# (xlabel, ylabel, etc. are usually retained)
plot 'density_Z_PRH.xvg' using 1:2 with lines linewidth 2 title 'Propanol'
unset output

set output 'density_profile_total.png'
set title "Total System Density Profile along Z-axis"
plot 'density_Z_total.xvg' using 1:2 with lines linewidth 2 title 'Total System'
unset output

# Overlaying MOH and PRH (most informative):
set output 'density_profile_overlay_MOH_PRH.png'
set title "Component Density Profiles along Z-axis"
plot 'density_Z_MOH.xvg' using 1:2 with lines linewidth 2 title 'Methanol (MOH)', \
     'density_Z_PRH.xvg' using 1:2 with lines linewidth 2 title 'Propanol (PRH)'
unset output

# Exit gnuplot:
# exit
\end{lstlisting}
Download these \texttt{.png} files to view them.

\subsubsection*{Interpreting the Density Profile Plots}
\begin{itemize}
    \item \textbf{Identify Phases:} Clearly see a high-density region (liquid slab) in the center and two low-density regions (vapor phases) at the ends.
    Example: If Z-box is 0-12 nm, liquid might be Z=3-9 nm, vapor Z=0-3 nm and Z=9-12 nm.
    \item \textbf{Check for Symmetry:} Ideally, two vapor phases (and interfaces) should be roughly symmetrical if well-equilibrated.
    \item \textbf{Flat Plateaus:} Density within bulk liquid and vapor regions should be relatively flat. A consistent slope might indicate incomplete equilibration or non-bulk region analysis.
    \item \textbf{Interface Width:} The transition region width is also a property of interest.
\end{itemize}

\begin{figure}[h!]
\centering
\includegraphics[width=0.7\textwidth]{density_profile_overlay.png} % Replace with actual filename
% \fbox{\parbox[c][6cm][c]{0.7\textwidth}{\centering Placeholder for Density Profile Overlay Image (your\_density\_profile\_overlay.pn)}}
\caption{Example density profiles for Methanol and Propanol along the Z-axis, clearly showing the higher density liquid slab in the center and the low-density vapor phases at the ends of the box. The Z-axis represents the dimension perpendicular to the liquid-vapor interfaces.}
\label{fig:density_profile_example}
\end{figure}

This visual inspection is crucial. If no clear phase separation, the simulation might need to run longer or other issues exist. Assuming clear phases, proceed to quantitative analysis.

\subsection{5.2. Calculating Phase Compositions \& Vapor Pressure}
\label{ssec:phase_comp_vp}
From density profiles and the energy file, calculate VLE properties: compositions of coexisting liquid/vapor phases, and vapor pressure.

\subsubsection*{1. Determining Equilibrium Densities and Phase Boundaries}
Examine density profile plots (e.g., \texttt{density\_Z\_total.xvg}):
\begin{itemize}
    \item Identify Z-coordinates for the bulk liquid phase (center, flat plateau). Example: $Z_{\text{liquid\_start}} = 4.0$~nm, $Z_{\text{liquid\_end}} = 8.0$~nm.
    \item Identify Z-coordinates for bulk vapor phases (ends, low density, flat plateaus). Example: $Z_{\text{vapor1\_start}} = 0.5$~nm, $Z_{\text{vapor1\_end}} = 2.5$~nm AND $Z_{\text{vapor2\_start}} = 9.5$~nm, $Z_{\text{vapor2\_end}} = 11.5$~nm. Average over both vapor regions if symmetrical.
    \item Choose regions away from rapidly changing interfaces.
\end{itemize}

\subsubsection*{Extract Average Densities from \texttt{.xvg} Files}
The \texttt{.xvg} files contain Z-coordinate (col 1) and density (kg/m$^3$, col 2). Calculate average density for each component (\texttt{MOH}, \texttt{PRH}) and total system within identified liquid/vapor regions.

\textbf{Method A: Manual/Spreadsheet (for a few points):}
Open \texttt{.xvg} (e.g., \texttt{density\_Z\_MOH.xvg}), identify lines for chosen Z-ranges, copy density values (col 2), and average.

\textbf{Method B: Scripting (Recommended):}
Use Python, Awk, etc., to read \texttt{.xvg}, filter lines by Z-range, and average densities.
Conceptual Python Example:
\begin{lstlisting}[language=python, caption=Python script to calculate average density from .xvg, label=lst:avg_density_py]
def get_average_density(xvg_filepath, z_start, z_end):
    densities_in_range = []
    with open(xvg_filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue # Skip comments
            parts = line.split()
            if len(parts) < 2: continue # Skip malformed lines
            try:
                z_coord = float(parts[0])
                density_val = float(parts[1]) # Corrected variable name
                if z_start <= z_coord <= z_end:
                    densities_in_range.append(density_val)
            except ValueError:
                continue # Skip lines where conversion to float fails
    if not densities_in_range:
        return 0.0 # Or raise an error, or return None
    return sum(densities_in_range) / len(densities_in_range)

# Define your Z-ranges based on visual inspection of plots
z_liq_start, z_liq_end = 4.0, 8.0
z_vap1_start, z_vap1_end = 0.5, 2.5
z_vap2_start, z_vap2_end = 9.5, 11.5

# Calculate average densities
avg_rho_MOH_liq = get_average_density('density_Z_MOH.xvg', z_liq_start, z_liq_end)
avg_rho_PRH_liq = get_average_density('density_Z_PRH.xvg', z_liq_start, z_liq_end)
# ... and for vapor phases (average of two vapor regions if used)
avg_rho_MOH_vap_1 = get_average_density('density_Z_MOH.xvg', z_vap1_start, z_vap1_end)
avg_rho_MOH_vap_2 = get_average_density('density_Z_MOH.xvg', z_vap2_start, z_vap2_end)
avg_rho_MOH_vap = (avg_rho_MOH_vap_1 + avg_rho_MOH_vap_2) / 2.0 if avg_rho_MOH_vap_1 is not None and avg_rho_MOH_vap_2 is not None else 0.0

avg_rho_PRH_vap_1 = get_average_density('density_Z_PRH.xvg', z_vap1_start, z_vap1_end)
avg_rho_PRH_vap_2 = get_average_density('density_Z_PRH.xvg', z_vap2_start, z_vap2_end)
avg_rho_PRH_vap = (avg_rho_PRH_vap_1 + avg_rho_PRH_vap_2) / 2.0 if avg_rho_PRH_vap_1 is not None and avg_rho_PRH_vap_2 is not None else 0.0


print(f"Avg Liquid MOH density: {avg_rho_MOH_liq:.4f} kg/m^3")
print(f"Avg Liquid PRH density: {avg_rho_PRH_liq:.4f} kg/m^3")
print(f"Avg Vapor MOH density: {avg_rho_MOH_vap:.4f} kg/m^3")
print(f"Avg Vapor PRH density: {avg_rho_PRH_vap:.4f} kg/m^3")
\end{lstlisting}
You will obtain: $\langle\rho_{\text{MOH,liquid}}\rangle$, $\langle\rho_{\text{PRH,liquid}}\rangle$, $\langle\rho_{\text{MOH,vapor}}\rangle$, $\langle\rho_{\text{PRH,vapor}}\rangle$.

\subsubsection*{2. Calculate Mole Fractions in Each Phase ($x_i, y_i$)}
Convert average mass densities to mole fractions.
\textbf{Molar Masses (M):}
\begin{itemize}
    \item Methanol (CH$_3$OH): $M_{\text{MOH}} \approx 32.04$~g/mol = $0.03204$~kg/mol
    \item 1-Propanol (CH$_3$CH$_2$CH$_2$OH): $M_{\text{PRH}} \approx 60.10$~g/mol = $0.06010$~kg/mol
\end{itemize}
\textbf{Calculate Molar Concentrations ($c = \rho / M$) in mol/m$^3$:}
\begin{align*}
    c_{\text{MOH,liquid}} &= \langle\rho_{\text{MOH,liquid}}\rangle / M_{\text{MOH}} \\
    c_{\text{PRH,liquid}} &= \langle\rho_{\text{PRH,liquid}}\rangle / M_{\text{PRH}} \\
    c_{\text{MOH,vapor}} &= \langle\rho_{\text{MOH,vapor}}\rangle / M_{\text{MOH}} \\
    c_{\text{PRH,vapor}} &= \langle\rho_{\text{PRH,vapor}}\rangle / M_{\text{PRH}}
\end{align*}
\textbf{Calculate Mole Fractions:}
\begin{itemize}
    \item Liquid Phase:
    \begin{align*}
        x_{\text{MOH}} &= c_{\text{MOH,liquid}} / (c_{\text{MOH,liquid}} + c_{\text{PRH,liquid}}) \\
        x_{\text{PRH}} &= c_{\text{PRH,liquid}} / (c_{\text{MOH,liquid}} + c_{\text{PRH,liquid}}) \\
        &\text{(Check: } x_{\text{MOH}} + x_{\text{PRH}} \approx 1)
    \end{align*}
    \item Vapor Phase:
    \begin{align*}
        y_{\text{MOH}} &= c_{\text{MOH,vapor}} / (c_{\text{MOH,vapor}} + c_{\text{PRH,vapor}}) \\
        y_{\text{PRH}} &= c_{\text{PRH,vapor}} / (c_{\text{MOH,vapor}} + c_{\text{PRH,vapor}}) \\
        &\text{(Check: } y_{\text{MOH}} + y_{\text{PRH}} \approx 1)
    \end{align*}
\end{itemize}
You now have equilibrium mole fractions ($x_{\text{MOH}}, y_{\text{MOH}}$) for methanol.

\subsubsection*{3. Determine the Vapor Pressure ($P_{\text{vap}}$)}
The VLE production run was NVT. Average pressure in the equilibrated vapor phase is $P_{\text{vap}}$.
Use \texttt{gmx energy}:
\begin{lstlisting}[language=bash, caption=Extracting pressure data with \texttt{gmx energy}, label=lst:gmx_energy_pressure]
gmx energy -f system_vle.edr -o pressure_vle.xvg
\end{lstlisting}
Select "Pressure" (and optionally "Time").
\textbf{Analyze \texttt{pressure\_vle.xvg}:}
Plot Pressure vs. Time. Discard initial non-equilibrated portion. Calculate average pressure over the stable, equilibrated portion. This is $P_{\text{vap}}$. GROMACS usually outputs pressure in bar.
Conceptual Python Example:
\begin{lstlisting}[language=python, caption=Python script to calculate average pressure from .xvg, label=lst:avg_pressure_py]
def get_average_pressure(xvg_filepath, time_start_ps):
    pressures_in_range = []
    with open(xvg_filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) < 2: continue
            try:
                time_ps = float(parts[0])
                pressure_val = float(parts[1]) # Corrected variable name
                if time_ps >= time_start_ps:
                    pressures_in_range.append(pressure_val)
            except ValueError:
                continue
    if not pressures_in_range:
        return 0.0 # Or raise an error, or return None
    return sum(pressures_in_range) / len(pressures_in_range)

# Example: Start averaging after 5 ns (5000 ps) of VLE run
avg_Pvap = get_average_pressure('pressure_vle.xvg', 5000.0)
print(f"Average Vapor Pressure (Pvap): {avg_Pvap:.4f} bar")
\end{lstlisting}

\subsubsection*{Summary of Data for One VLE Point}
From one VLE simulation at temperature $T$, you have:
\begin{itemize}
    \item Liquid phase mole fraction of methanol: $x_{\text{MOH}}$
    \item Vapor phase mole fraction of methanol: $y_{\text{MOH}}$
    \item Vapor pressure of the mixture: $P_{\text{vap}}$
\end{itemize}
This set ($T, P_{\text{vap}}, x_{\text{MOH}}, y_{\text{MOH}}$) is one tie-line.

\subsubsection*{Uncertainty Estimation (Advanced)}
For rigorous work, estimate uncertainties:
\begin{itemize}
    \item \textbf{Block Averaging:} Divide equilibrated trajectory into blocks, calculate property for each block, then standard error of the mean.
    \item Analyzing fluctuations over time.
\end{itemize}
This provides error bars for VLE data points.

\subsection{5.3. Building VLE Diagrams}
\label{ssec:vle_diagrams}
Phases \ref{sec:phase1}-\ref{sec:phase4} and Sections \ref{ssec:gmx_density_profiles}-\ref{ssec:phase_comp_vp} detail a single VLE simulation at one temperature and initial composition, yielding ($T, P_{\text{vap}}, x_i, y_i$). To construct VLE diagrams (T-xy, P-xy, xy), repeat the process, varying temperature or initial composition.

\subsubsection*{1. Generating Data for a T-xy Diagram}
Shows bubble/dew point temperatures vs. composition at fixed external pressure.
\textbf{Strategy with NVT Slab Simulations:}
NVT slab simulations determine their own $P_{\text{vap}}$. To approximate a T-xy diagram at a target pressure (e.g., 1 atm = 1.01325 bar):
\begin{enumerate}
    \item Choose temperatures ($T_1, T_2, \dots$) spanning the expected boiling range.
    \item For each $T_i$:
    \begin{itemize}
        \item Perform full VLE workflow (equilibration, slab, NVT VLE production at $T_i$, analysis).
        \item Obtain ($x_{\text{MOH},i}, y_{\text{MOH},i}, P_{\text{vap},i}$).
    \end{itemize}
\end{enumerate}
\textbf{Data Collection Table:}
\begin{table}[h!]
\centering
\caption{Example data for T-xy diagram construction.}
\label{tab:txy_data}
\begin{tabular}{@{}cccc@{}}
\toprule
Temperature (K) ($T_i$) & $x_{\text{MOH},i}$ & $y_{\text{MOH},i}$ & $P_{\text{vap},i}$ (bar) \\ \midrule
$T_1$ & $x_1$ & $y_1$ & $P_1$ \\
$T_2$ & $x_2$ & $y_2$ & $P_2$ \\
$\dots$ & $\dots$ & $\dots$ & $\dots$ \\ \bottomrule
\end{tabular}
\end{table}
\textbf{Plotting T-xy Diagram:}
\begin{itemize}
    \item For exactly 1 atm: Find, for several overall compositions, the temperature where $P_{\text{vap}} = 1$~atm (complex, involves multiple simulations and interpolation).
    \item More direct plot from NVT slabs: Plot $T_i$ (y-axis) vs. $x_{\text{MOH},i}$ and $y_{\text{MOH},i}$ (x-axis).
    \begin{itemize}
        \item Bubble point curve: $T_i$ vs. $x_{\text{MOH},i}$.
        \item Dew point curve: $T_i$ vs. $y_{\text{MOH},i}$.
    \end{itemize}
    These curves represent VLE at self-generated $P_{\text{vap},i}$.
\end{itemize}

\subsubsection*{2. Generating Data for a P-xy Diagram (at Constant Temperature)}
Shows bubble/dew point pressures vs. composition at fixed temperature.
\textbf{Strategy:}
\begin{enumerate}
    \item Fix simulation temperature ($T$).
    \item Vary Overall Initial Composition (in \texttt{packmol.inp}, Section \ref{ssec:packmol_packing}). Examples:
    \begin{itemize}
        \item Run 1: High Methanol (e.g., 80\% MOH)
        \item Run 2: Medium Methanol (e.g., 50\% MOH)
        \item Run 3: Low Methanol (e.g., 20\% MOH)
    \end{itemize}
    \item For each initial composition:
    \begin{itemize}
        \item Perform full VLE workflow at fixed $T$.
        \item Obtain coexisting ($x_{\text{MOH},i}, y_{\text{MOH},i}$) and $P_{\text{vap},i}$.
    \end{itemize}
\end{enumerate}
\textbf{Data Collection Table (fixed $T$):}
\begin{table}[h!]
\centering
\caption{Example data for P-xy diagram construction (fixed T).}
\label{tab:pxy_data}
\begin{tabular}{@{}cccc@{}}
\toprule
Overall Start (approx.) & $x_{\text{MOH},i}$ (Equil. Liq.) & $y_{\text{MOH},i}$ (Equil. Vap.) & $P_{\text{vap},i}$ (bar) \\ \midrule
80\% MOH & $x_1$ & $y_1$ & $P_1$ \\
50\% MOH & $x_2$ & $y_2$ & $P_2$ \\
20\% MOH & $x_3$ & $y_3$ & $P_3$ \\
$\dots$ & $\dots$ & $\dots$ & $\dots$ \\ \bottomrule
\end{tabular}
\end{table}
\textbf{Plotting P-xy Diagram:} Plot $P_{\text{vap},i}$ (y-axis) vs. $x_{\text{MOH},i}$ and $y_{\text{MOH},i}$ (x-axis).

\subsubsection*{3. Generating Data for an xy-Diagram (Relative Volatility Plot)}
Plots vapor mole fraction ($y_i$) vs. liquid mole fraction ($x_i$) at constant T or P. Useful for separability.
\textbf{Strategy:}
\begin{itemize}
    \item Use data from T-xy or P-xy generation.
    \item For each tie-line, you have ($x_{\text{MOH}}, y_{\text{MOH}}$).
    \item Plot $y_{\text{MOH}}$ (y-axis) vs. $x_{\text{MOH}}$ (x-axis).
    \item Often, a diagonal line ($y=x$) is drawn. Curve above diagonal $\implies$ component is more volatile.
\end{itemize}

\subsubsection*{Visualizing Diagrams (e.g., Gnuplot)}
Conceptual Gnuplot for T-xy (data in \texttt{vle\_data.dat}: Temp $x_{\text{MOH}}$ $y_{\text{MOH}}$ $P_{\text{vap}}$):
\begin{lstlisting}[language=gnuplot, caption=Conceptual Gnuplot script for T-xy diagram, label=lst:gnuplot_txy]
set title "T-xy Diagram for Methanol-Propanol"
set xlabel "Mole Fraction Methanol (x_{MOH}, y_{MOH})"
set ylabel "Temperature (K)"
set xrange [0:1]
# set yrange [min_temp:max_temp] # Adjust to your temperature range
set grid
set key top right

plot 'vle_data.dat' using 2:1 with linespoints title 'Bubble Point (T vs x_{MOH})', \
     'vle_data.dat' using 3:1 with linespoints title 'Dew Point (T vs y_{MOH})'
\end{lstlisting}

\subsubsection*{Key Considerations for Building Diagrams}
\begin{itemize}
    \item \textbf{Number of Points:} 5-7 points across the range is a good start.
    \item \textbf{Equilibration:} Ensure each VLE simulation is truly equilibrated (stable density profiles, compositions, pressure).
    \item \textbf{Computational Cost:} Each point requires a full, lengthy simulation.
    \item \textbf{Force Field Limitations:} Accuracy depends on force field quality (GAFF2 here). Compare with experimental data if available.
\end{itemize}

\subsubsection*{Your Next Steps (Practical)}
\begin{enumerate}
    \item Complete analysis for current 298.15 K simulation:
    \begin{itemize}
        \item Determine Z-ranges for bulk liquid/vapor.
        \item Extract average densities for \texttt{MOH} and \texttt{PRH}.
        \item Calculate $x_{\text{MOH}}$ and $y_{\text{MOH}}$.
        \item Calculate $P_{\text{vap}}$ from \texttt{system\_vle.edr}.
    \end{itemize}
    You have one VLE data point: ($T=298.15$~K, $P_{\text{vap}}, x_{\text{MOH}}, y_{\text{MOH}}$).
    \item Plan Your Next Simulation:
    \begin{itemize}
        \item If T-xy: Choose new temperature (e.g., 310 K or 320 K). Update \texttt{ref\_t} in \texttt{nvt\_equil.mdp}, \texttt{npt.mdp}, and \texttt{vle\_nvt.mdp}.
        \item If P-xy: Keep $T=298.15$~K. Change overall composition in \texttt{packmol.inp} (molecule numbers).
    \end{itemize}
\end{enumerate}
This systematic approach builds VLE diagrams.

% --- END OF PHASE 4 --- (Assumes previous content for Phases 1, 2, 3, and 4 is in place)

%==================================================================
%   PHASE 5: ANALYSIS AND INTERPRETATION
%==================================================================
\section[Phase 5: Analysis and Interpretation]{Phase 5: Analysis and Interpretation} % Manual line break for ToC
\label{sec:phase5}
After your long VLE production run (\texttt{gmx mdrun -deffnm system\_vle}) has completed successfully, you will have a set of output files (primarily \texttt{system\_vle.xtc}, \texttt{system\_vle.edr}, and \texttt{system\_vle.tpr}) that contain the raw data about your system's behavior at the simulated temperature. The goal of this phase is to analyze this data to determine the properties of the coexisting liquid and vapor phases, which constitute one point on your VLE diagram.

\subsection{5.1. Density Profiles (\texttt{gmx density})}
\label{ssec:gmx_density_profiles}

\subsubsection*{Purpose}
The most direct way to confirm the formation of distinct liquid and vapor phases and to determine their boundaries and compositions is by calculating the density profile of each component (and the total system) along the Z-axis (the dimension perpendicular to the liquid-vapor interfaces).
\begin{itemize}
    \item A high-density plateau in the center of the box indicates the liquid slab.
    \item Low-density plateaus at the ends of the box indicate the vapor phases.
    \item The regions of transition between these plateaus are the liquid-vapor interfaces.
\end{itemize}

\subsubsection*{The \texttt{gmx density} Tool}
GROMACS provides the \texttt{gmx density} tool for this purpose.

\subsubsection*{Commands to Generate Density Profiles}
You will run \texttt{gmx density} separately for each component and for the total system. It's important to analyze a significant portion of your trajectory, potentially discarding an initial part if the system was still evolving towards phase equilibrium (e.g., if interface formation took some nanoseconds). The \texttt{-b} (begin time) and \texttt{-e} (end time) flags in \texttt{gmx density} can be used for this, but for a first pass, analyzing the whole trajectory is common.

\textbf{For Methanol (\texttt{MOH}):}
\begin{lstlisting}[language=bash, caption=Generating methanol density profile with \texttt{gmx density}, label=lst:gmx_density_moh]
gmx density -f system_vle.xtc -s system_vle.tpr -d Z -o density_Z_MOH.xvg -sl 120
\end{lstlisting}
\begin{itemize}
    \item \texttt{-f system\_vle.xtc}: Input trajectory file.
    \item \texttt{-s system\_vle.tpr}: Input run file (topology and system information).
    \item \texttt{-d Z}: Calculate density along the Z-axis.
    \item \texttt{-o density\_Z\_MOH.xvg}: Output file for methanol density profile (\texttt{.xvg} is Grace/XMGRACE format).
    \item \texttt{-sl 120}: Number of slices along Z-axis. For a $\approx 12$~nm box, 120 slices give 0.1~nm resolution. Adjust for reasonable resolution (0.05 to 0.1~nm/slice often good).
\end{itemize}
\textbf{Interactive Prompt:} \texttt{gmx density} will list atom groups (e.g., System, Other, MOH, PRH). Select the group for Methanol (e.g., by typing its number). In our case, \texttt{MOH} was group 2.

\textbf{For Propanol (\texttt{PRH}):}
\begin{lstlisting}[language=bash, caption=Generating propanol density profile with \texttt{gmx density}, label=lst:gmx_density_prh]
gmx density -f system_vle.xtc -s system_vle.tpr -d Z -o density_Z_PRH.xvg -sl 120
\end{lstlisting}
When prompted, select the group for Propanol (e.g., \texttt{PRH}, group 3).

\textbf{For the Total System:}
\begin{lstlisting}[language=bash, caption=Generating total system density profile with \texttt{gmx density}, label=lst:gmx_density_total]
gmx density -f system_vle.xtc -s system_vle.tpr -d Z -o density_Z_total.xvg -sl 120
\end{lstlisting}
When prompted, select the group for the entire System (usually group 0 or 1).

\subsubsection*{Output \texttt{.xvg} Files}
Each command produces an \texttt{.xvg} file (plain text):
\begin{itemize}
    \item Lines starting with \texttt{\#} or \texttt{@} are comments/metadata.
    \item Data lines:
    \begin{itemize}
        \item Column 1: Z-coordinate (nm) of the slice center.
        \item Column 2: Density (kg/m$^3$) in that slice, averaged over analyzed frames.
    \end{itemize}
\end{itemize}

\subsubsection*{Visualizing Density Profiles (e.g., with Gnuplot)}
Use Gnuplot to visualize profiles and save images.
\begin{lstlisting}[language=gnuplot, caption=Gnuplot script for visualizing density profiles, label=lst:gnuplot_density]
# Start gnuplot:
# gnuplot

# Plotting individual profiles and saving:
set terminal pngcairo size 800,600 enhanced font "arial,10"
set output 'density_profile_MOH.png'
set title "Density Profile of Methanol (MOH) along Z-axis"
set xlabel "Z-coordinate (nm)"
set ylabel "Density (kg/m^3)"
set grid
set yrange [0:*] # Ensure y-axis starts at 0, auto-scale max
plot 'density_Z_MOH.xvg' using 1:2 with lines linewidth 2 title 'Methanol'
unset output

set output 'density_profile_PRH.png'
set title "Density Profile of Propanol (PRH) along Z-axis"
# (xlabel, ylabel, etc. are usually retained)
plot 'density_Z_PRH.xvg' using 1:2 with lines linewidth 2 title 'Propanol'
unset output

set output 'density_profile_total.png'
set title "Total System Density Profile along Z-axis"
plot 'density_Z_total.xvg' using 1:2 with lines linewidth 2 title 'Total System'
unset output

# Overlaying MOH and PRH (most informative):
set output 'density_profile_overlay_MOH_PRH.png'
set title "Component Density Profiles along Z-axis"
plot 'density_Z_MOH.xvg' using 1:2 with lines linewidth 2 title 'Methanol (MOH)', \
     'density_Z_PRH.xvg' using 1:2 with lines linewidth 2 title 'Propanol (PRH)'
unset output

# Exit gnuplot:
# exit
\end{lstlisting}
Download these \texttt{.png} files to view them.

\subsubsection*{Interpreting the Density Profile Plots}
\begin{itemize}
    \item \textbf{Identify Phases:} Clearly see a high-density region (liquid slab) in the center and two low-density regions (vapor phases) at the ends.
    Example: If Z-box is 0-12 nm, liquid might be Z=3-9 nm, vapor Z=0-3 nm and Z=9-12 nm.
    \item \textbf{Check for Symmetry:} Ideally, two vapor phases (and interfaces) should be roughly symmetrical if well-equilibrated.
    \item \textbf{Flat Plateaus:} Density within bulk liquid and vapor regions should be relatively flat. A consistent slope might indicate incomplete equilibration or non-bulk region analysis.
    \item \textbf{Interface Width:} The transition region width is also a property of interest.
\end{itemize}

\begin{figure}[h!]
\centering
% \includegraphics[width=0.7\textwidth]{your_density_profile_overlay.png} % Replace with actual filename
\fbox{\parbox[c][6cm][c]{0.7\textwidth}{\centering Placeholder for Density Profile Overlay Image (your\_density\_profile\_overlay.png)}}
\caption{Example density profiles for Methanol and Propanol along the Z-axis, clearly showing the higher density liquid slab in the center and the low-density vapor phases at the ends of the box. The Z-axis represents the dimension perpendicular to the liquid-vapor interfaces.}
\label{fig:density_profile_example}
\end{figure}

This visual inspection is crucial. If no clear phase separation, the simulation might need to run longer or other issues exist. Assuming clear phases, proceed to quantitative analysis.

\subsection{5.2. Calculating Phase Compositions \& Vapor Pressure}
\label{ssec:phase_comp_vp}
From density profiles and the energy file, calculate VLE properties: compositions of coexisting liquid/vapor phases, and vapor pressure.

\subsubsection*{1. Determining Equilibrium Densities and Phase Boundaries}
Examine density profile plots (e.g., \texttt{density\_Z\_total.xvg}):
\begin{itemize}
    \item Identify Z-coordinates for the bulk liquid phase (center, flat plateau). Example: $Z_{\text{liquid\_start}} = 4.0$~nm, $Z_{\text{liquid\_end}} = 8.0$~nm.
    \item Identify Z-coordinates for bulk vapor phases (ends, low density, flat plateaus). Example: $Z_{\text{vapor1\_start}} = 0.5$~nm, $Z_{\text{vapor1\_end}} = 2.5$~nm AND $Z_{\text{vapor2\_start}} = 9.5$~nm, $Z_{\text{vapor2\_end}} = 11.5$~nm. Average over both vapor regions if symmetrical.
    \item Choose regions away from rapidly changing interfaces.
\end{itemize}

\subsubsection*{Extract Average Densities from \texttt{.xvg} Files}
The \texttt{.xvg} files contain Z-coordinate (col 1) and density (kg/m$^3$, col 2). Calculate average density for each component (\texttt{MOH}, \texttt{PRH}) and total system within identified liquid/vapor regions.

\textbf{Method A: Manual/Spreadsheet (for a few points):}
Open \texttt{.xvg} (e.g., \texttt{density\_Z\_MOH.xvg}), identify lines for chosen Z-ranges, copy density values (col 2), and average.

\textbf{Method B: Scripting (Recommended):}
Use Python, Awk, etc., to read \texttt{.xvg}, filter lines by Z-range, and average densities.
Conceptual Python Example:
\begin{lstlisting}[language=python, caption=Python script to calculate average density from .xvg, label=lst:avg_density_py]
def get_average_density(xvg_filepath, z_start, z_end):
    densities_in_range = []
    with open(xvg_filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue # Skip comments
            parts = line.split()
            if len(parts) < 2: continue # Skip malformed lines
            try:
                z_coord = float(parts[0])
                density_val = float(parts[1]) # Corrected variable name
                if z_start <= z_coord <= z_end:
                    densities_in_range.append(density_val)
            except ValueError:
                continue # Skip lines where conversion to float fails
    if not densities_in_range:
        return 0.0 # Or raise an error, or return None
    return sum(densities_in_range) / len(densities_in_range)

# Define your Z-ranges based on visual inspection of plots
z_liq_start, z_liq_end = 4.0, 8.0
z_vap1_start, z_vap1_end = 0.5, 2.5
z_vap2_start, z_vap2_end = 9.5, 11.5

# Calculate average densities
avg_rho_MOH_liq = get_average_density('density_Z_MOH.xvg', z_liq_start, z_liq_end)
avg_rho_PRH_liq = get_average_density('density_Z_PRH.xvg', z_liq_start, z_liq_end)
# ... and for vapor phases (average of two vapor regions if used)
avg_rho_MOH_vap_1 = get_average_density('density_Z_MOH.xvg', z_vap1_start, z_vap1_end)
avg_rho_MOH_vap_2 = get_average_density('density_Z_MOH.xvg', z_vap2_start, z_vap2_end)
avg_rho_MOH_vap = (avg_rho_MOH_vap_1 + avg_rho_MOH_vap_2) / 2.0 if avg_rho_MOH_vap_1 is not None and avg_rho_MOH_vap_2 is not None else 0.0

avg_rho_PRH_vap_1 = get_average_density('density_Z_PRH.xvg', z_vap1_start, z_vap1_end)
avg_rho_PRH_vap_2 = get_average_density('density_Z_PRH.xvg', z_vap2_start, z_vap2_end)
avg_rho_PRH_vap = (avg_rho_PRH_vap_1 + avg_rho_PRH_vap_2) / 2.0 if avg_rho_PRH_vap_1 is not None and avg_rho_PRH_vap_2 is not None else 0.0


print(f"Avg Liquid MOH density: {avg_rho_MOH_liq:.4f} kg/m^3")
print(f"Avg Liquid PRH density: {avg_rho_PRH_liq:.4f} kg/m^3")
print(f"Avg Vapor MOH density: {avg_rho_MOH_vap:.4f} kg/m^3")
print(f"Avg Vapor PRH density: {avg_rho_PRH_vap:.4f} kg/m^3")
\end{lstlisting}
You will obtain: $\langle\rho_{\text{MOH,liquid}}\rangle$, $\langle\rho_{\text{PRH,liquid}}\rangle$, $\langle\rho_{\text{MOH,vapor}}\rangle$, $\langle\rho_{\text{PRH,vapor}}\rangle$.

\subsubsection*{2. Calculate Mole Fractions in Each Phase ($x_i, y_i$)}
Convert average mass densities to mole fractions.
\textbf{Molar Masses (M):}
\begin{itemize}
    \item Methanol (CH$_3$OH): $M_{\text{MOH}} \approx 32.04$~g/mol = $0.03204$~kg/mol
    \item 1-Propanol (CH$_3$CH$_2$CH$_2$OH): $M_{\text{PRH}} \approx 60.10$~g/mol = $0.06010$~kg/mol
\end{itemize}
\textbf{Calculate Molar Concentrations ($c = \rho / M$) in mol/m$^3$:}
\begin{align*}
    c_{\text{MOH,liquid}} &= \langle\rho_{\text{MOH,liquid}}\rangle / M_{\text{MOH}} \\
    c_{\text{PRH,liquid}} &= \langle\rho_{\text{PRH,liquid}}\rangle / M_{\text{PRH}} \\
    c_{\text{MOH,vapor}} &= \langle\rho_{\text{MOH,vapor}}\rangle / M_{\text{MOH}} \\
    c_{\text{PRH,vapor}} &= \langle\rho_{\text{PRH,vapor}}\rangle / M_{\text{PRH}}
\end{align*}
\textbf{Calculate Mole Fractions:}
\begin{itemize}
    \item Liquid Phase:
    \begin{align*}
        x_{\text{MOH}} &= c_{\text{MOH,liquid}} / (c_{\text{MOH,liquid}} + c_{\text{PRH,liquid}}) \\
        x_{\text{PRH}} &= c_{\text{PRH,liquid}} / (c_{\text{MOH,liquid}} + c_{\text{PRH,liquid}}) \\
        &\text{(Check: } x_{\text{MOH}} + x_{\text{PRH}} \approx 1)
    \end{align*}
    \item Vapor Phase:
    \begin{align*}
        y_{\text{MOH}} &= c_{\text{MOH,vapor}} / (c_{\text{MOH,vapor}} + c_{\text{PRH,vapor}}) \\
        y_{\text{PRH}} &= c_{\text{PRH,vapor}} / (c_{\text{MOH,vapor}} + c_{\text{PRH,vapor}}) \\
        &\text{(Check: } y_{\text{MOH}} + y_{\text{PRH}} \approx 1)
    \end{align*}
\end{itemize}
You now have equilibrium mole fractions ($x_{\text{MOH}}, y_{\text{MOH}}$) for methanol.

\subsubsection*{3. Determine the Vapor Pressure ($P_{\text{vap}}$)}
The VLE production run was NVT. Average pressure in the equilibrated vapor phase is $P_{\text{vap}}$.
Use \texttt{gmx energy}:
\begin{lstlisting}[language=bash, caption=Extracting pressure data with \texttt{gmx energy}, label=lst:gmx_energy_pressure]
gmx energy -f system_vle.edr -o pressure_vle.xvg
\end{lstlisting}
Select "Pressure" (and optionally "Time").
\textbf{Analyze \texttt{pressure\_vle.xvg}:}
Plot Pressure vs. Time. Discard initial non-equilibrated portion. Calculate average pressure over the stable, equilibrated portion. This is $P_{\text{vap}}$. GROMACS usually outputs pressure in bar.
Conceptual Python Example:
\begin{lstlisting}[language=python, caption=Python script to calculate average pressure from .xvg, label=lst:avg_pressure_py]
def get_average_pressure(xvg_filepath, time_start_ps):
    pressures_in_range = []
    with open(xvg_filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) < 2: continue
            try:
                time_ps = float(parts[0])
                pressure_val = float(parts[1]) # Corrected variable name
                if time_ps >= time_start_ps:
                    pressures_in_range.append(pressure_val)
            except ValueError:
                continue
    if not pressures_in_range:
        return 0.0 # Or raise an error, or return None
    return sum(pressures_in_range) / len(pressures_in_range)

# Example: Start averaging after 5 ns (5000 ps) of VLE run
avg_Pvap = get_average_pressure('pressure_vle.xvg', 5000.0)
print(f"Average Vapor Pressure (Pvap): {avg_Pvap:.4f} bar")
\end{lstlisting}

\subsubsection*{Summary of Data for One VLE Point}
From one VLE simulation at temperature $T$, you have:
\begin{itemize}
    \item Liquid phase mole fraction of methanol: $x_{\text{MOH}}$
    \item Vapor phase mole fraction of methanol: $y_{\text{MOH}}$
    \item Vapor pressure of the mixture: $P_{\text{vap}}$
\end{itemize}
This set ($T, P_{\text{vap}}, x_{\text{MOH}}, y_{\text{MOH}}$) is one tie-line.

\subsubsection*{Uncertainty Estimation (Advanced)}
For rigorous work, estimate uncertainties:
\begin{itemize}
    \item \textbf{Block Averaging:} Divide equilibrated trajectory into blocks, calculate property for each block, then standard error of the mean.
    \item Analyzing fluctuations over time.
\end{itemize}
This provides error bars for VLE data points.

\subsection{5.3. Building VLE Diagrams}
\label{ssec:vle_diagrams}
Phases \ref{sec:phase1}-\ref{sec:phase4} and Sections \ref{ssec:gmx_density_profiles}-\ref{ssec:phase_comp_vp} detail a single VLE simulation at one temperature and initial composition, yielding ($T, P_{\text{vap}}, x_i, y_i$). To construct VLE diagrams (T-xy, P-xy, xy), repeat the process, varying temperature or initial composition.

\subsubsection*{1. Generating Data for a T-xy Diagram}
Shows bubble/dew point temperatures vs. composition at fixed external pressure.
\textbf{Strategy with NVT Slab Simulations:}
NVT slab simulations determine their own $P_{\text{vap}}$. To approximate a T-xy diagram at a target pressure (e.g., 1 atm = 1.01325 bar):
\begin{enumerate}
    \item Choose temperatures ($T_1, T_2, \dots$) spanning the expected boiling range.
    \item For each $T_i$:
    \begin{itemize}
        \item Perform full VLE workflow (equilibration, slab, NVT VLE production at $T_i$, analysis).
        \item Obtain ($x_{\text{MOH},i}, y_{\text{MOH},i}, P_{\text{vap},i}$).
    \end{itemize}
\end{enumerate}
\textbf{Data Collection Table:}
\begin{table}[h!]
\centering
\caption{Example data for T-xy diagram construction.}
\label{tab:txy_data}
\begin{tabular}{@{}cccc@{}}
\toprule
Temperature (K) ($T_i$) & $x_{\text{MOH},i}$ & $y_{\text{MOH},i}$ & $P_{\text{vap},i}$ (bar) \\ \midrule
$T_1$ & $x_1$ & $y_1$ & $P_1$ \\
$T_2$ & $x_2$ & $y_2$ & $P_2$ \\
$\dots$ & $\dots$ & $\dots$ & $\dots$ \\ \bottomrule
\end{tabular}
\end{table}
\textbf{Plotting T-xy Diagram:}
\begin{itemize}
    \item For exactly 1 atm: Find, for several overall compositions, the temperature where $P_{\text{vap}} = 1$~atm (complex, involves multiple simulations and interpolation).
    \item More direct plot from NVT slabs: Plot $T_i$ (y-axis) vs. $x_{\text{MOH},i}$ and $y_{\text{MOH},i}$ (x-axis).
    \begin{itemize}
        \item Bubble point curve: $T_i$ vs. $x_{\text{MOH},i}$.
        \item Dew point curve: $T_i$ vs. $y_{\text{MOH},i}$.
    \end{itemize}
    These curves represent VLE at self-generated $P_{\text{vap},i}$.
\end{itemize}

\subsubsection*{2. Generating Data for a P-xy Diagram (at Constant Temperature)}
Shows bubble/dew point pressures vs. composition at fixed temperature.
\textbf{Strategy:}
\begin{enumerate}
    \item Fix simulation temperature ($T$).
    \item Vary Overall Initial Composition (in \texttt{packmol.inp}, Section \ref{ssec:packmol_packing}). Examples:
    \begin{itemize}
        \item Run 1: High Methanol (e.g., 80\% MOH)
        \item Run 2: Medium Methanol (e.g., 50\% MOH)
        \item Run 3: Low Methanol (e.g., 20\% MOH)
    \end{itemize}
    \item For each initial composition:
    \begin{itemize}
        \item Perform full VLE workflow at fixed $T$.
        \item Obtain coexisting ($x_{\text{MOH},i}, y_{\text{MOH},i}$) and $P_{\text{vap},i}$.
    \end{itemize}
\end{enumerate}
\textbf{Data Collection Table (fixed $T$):}
\begin{table}[h!]
\centering
\caption{Example data for P-xy diagram construction (fixed T).}
\label{tab:pxy_data}
\begin{tabular}{@{}cccc@{}}
\toprule
Overall Start (approx.) & $x_{\text{MOH},i}$ (Equil. Liq.) & $y_{\text{MOH},i}$ (Equil. Vap.) & $P_{\text{vap},i}$ (bar) \\ \midrule
80\% MOH & $x_1$ & $y_1$ & $P_1$ \\
50\% MOH & $x_2$ & $y_2$ & $P_2$ \\
20\% MOH & $x_3$ & $y_3$ & $P_3$ \\
$\dots$ & $\dots$ & $\dots$ & $\dots$ \\ \bottomrule
\end{tabular}
\end{table}
\textbf{Plotting P-xy Diagram:} Plot $P_{\text{vap},i}$ (y-axis) vs. $x_{\text{MOH},i}$ and $y_{\text{MOH},i}$ (x-axis).

\subsubsection*{3. Generating Data for an xy-Diagram (Relative Volatility Plot)}
Plots vapor mole fraction ($y_i$) vs. liquid mole fraction ($x_i$) at constant T or P. Useful for separability.
\textbf{Strategy:}
\begin{itemize}
    \item Use data from T-xy or P-xy generation.
    \item For each tie-line, you have ($x_{\text{MOH}}, y_{\text{MOH}}$).
    \item Plot $y_{\text{MOH}}$ (y-axis) vs. $x_{\text{MOH}}$ (x-axis).
    \item Often, a diagonal line ($y=x$) is drawn. Curve above diagonal $\implies$ component is more volatile.
\end{itemize}

\subsubsection*{Visualizing Diagrams (e.g., Gnuplot)}
Conceptual Gnuplot for T-xy (data in \texttt{vle\_data.dat}: Temp $x_{\text{MOH}}$ $y_{\text{MOH}}$ $P_{\text{vap}}$):
\begin{lstlisting}[language=gnuplot, caption=Conceptual Gnuplot script for T-xy diagram, label=lst:gnuplot_txy]
set title "T-xy Diagram for Methanol-Propanol"
set xlabel "Mole Fraction Methanol (x_{MOH}, y_{MOH})"
set ylabel "Temperature (K)"
set xrange [0:1]
# set yrange [min_temp:max_temp] # Adjust to your temperature range
set grid
set key top right

plot 'vle_data.dat' using 2:1 with linespoints title 'Bubble Point (T vs x_{MOH})', \
     'vle_data.dat' using 3:1 with linespoints title 'Dew Point (T vs y_{MOH})'
\end{lstlisting}

\subsubsection*{Key Considerations for Building Diagrams}
\begin{itemize}
    \item \textbf{Number of Points:} 5-7 points across the range is a good start.
    \item \textbf{Equilibration:} Ensure each VLE simulation is truly equilibrated (stable density profiles, compositions, pressure).
    \item \textbf{Computational Cost:} Each point requires a full, lengthy simulation.
    \item \textbf{Force Field Limitations:} Accuracy depends on force field quality (GAFF2 here). Compare with experimental data if available.
\end{itemize}

\subsubsection*{Your Next Steps (Practical)}
\begin{enumerate}
    \item Complete analysis for current 298.15 K simulation:
    \begin{itemize}
        \item Determine Z-ranges for bulk liquid/vapor.
        \item Extract average densities for \texttt{MOH} and \texttt{PRH}.
        \item Calculate $x_{\text{MOH}}$ and $y_{\text{MOH}}$.
        \item Calculate $P_{\text{vap}}$ from \texttt{system\_vle.edr}.
    \end{itemize}
    You have one VLE data point: ($T=298.15$~K, $P_{\text{vap}}, x_{\text{MOH}}, y_{\text{MOH}}$).
    \item Plan Your Next Simulation:
    \begin{itemize}
        \item If T-xy: Choose new temperature (e.g., 310 K or 320 K). Update \texttt{ref\_t} in \texttt{nvt\_equil.mdp}, \texttt{npt.mdp}, and \texttt{vle\_nvt.mdp}.
        \item If P-xy: Keep $T=298.15$~K. Change overall composition in \texttt{packmol.inp} (molecule numbers).
    \end{itemize}
\end{enumerate}
This systematic approach builds VLE diagrams.

% --- END OF PHASE 5 --- (Assumes previous content for Phases 1-5 is in place)

%==================================================================
%   TROUBLESHOOTING & COMMON ISSUES
%==================================================================
\section{Troubleshooting \& Common Issues}
\label{sec:troubleshooting}
This section is invaluable for anyone (including your future self) who encounters problems. It should be a "first aid" guide. For each issue, we describe the symptom, likely cause(s), and solution(s)/debugging steps.

\begin{itemize}
    \item \textbf{Issue: \texttt{antechamber} PDB Formatting Errors}
        \begin{itemize}
            \item \textit{Symptom:} \texttt{antechamber} exits with a fatal error related to PDB coordinate formatting.
            \item \textit{Likely Cause(s):} The input \texttt{.pdb} file for \texttt{antechamber} does not strictly adhere to PDB format, especially the fixed-column layout for atomic coordinates (X in 31-38, Y in 39-46, Z in 47-54, often F8.3). Spaces and justification matter.
            \item \textit{Solution(s) / Debugging Steps:}
            \begin{itemize}
                \item Manually inspect the input \texttt{.pdb} file with a text editor, paying close attention to column alignment for all \texttt{ATOM} or \texttt{HETATM} records.
                \item Use a text editor with a column mode or a ruler guide if possible.
                \item Ensure atom serial numbers, atom names, residue names, and residue sequence numbers are also correctly formatted and justified as per PDB spec, as misalignment here can shift coordinates.
                \item If the PDB was generated by a program, check its export options or try re-generating it. Using a reliable molecular editor (like Avogadro) to save the PDB often ensures correct formatting.
            \end{itemize}
        \end{itemize}

    \item \textbf{Issue: \texttt{antechamber} Internal Program Errors (\texttt{bondtype}, \texttt{atomtype})}
        \begin{itemize}
            \item \textit{Symptom:} \texttt{antechamber} reports a fatal error from an internal program like \texttt{bondtype} or \texttt{atomtype}. Subsequently, \texttt{parmchk2} (if run in the same script line) fails because the \texttt{.mol2} file was never created by \texttt{antechamber}.
            \item \textit{Likely Cause(s):} Often due to problematic input geometry in the PDB file (e.g., atoms too close, incorrect valencies implied by the structure) that \texttt{bondtype} cannot resolve.
            \item \textit{Solution(s) / Debugging Steps:}
            \begin{itemize}
                \item Carefully review the input \texttt{.pdb} file in a molecular visualizer.
                \item Perform a quick geometry optimization/cleanup in a tool like Avogadro before feeding the PDB to \texttt{antechamber}.
                \item Ensure all hydrogens are present and the molecule is chemically sensible.
            \end{itemize}
        \end{itemize}

    \item \textbf{Issue: \texttt{packmol} MOL2 Input Problems}
        \begin{itemize}
            \item \textit{Symptom:} \texttt{packmol} errors out when trying to read an input \texttt{.mol2} file, even if \texttt{filetype mol2} is specified in \texttt{packmol.inp}.
            \item \textit{Likely Cause(s):} This specific issue was encountered with Packmol v20.15.1. It seemed to ignore or misinterpret the \texttt{filetype mol2} directive and defaulted to trying to read the file as a PDB, which fails due to format differences.
            \item \textit{Solution(s) / Debugging Steps (Workaround Used):}
            \begin{itemize}
                \item Ensure initial PDB files (fed to \texttt{antechamber}) have globally unique atom names for each molecule type (e.g., \texttt{C1M, H1M} for methanol; \texttt{C1P, H2P} for propanol).
                \item In \texttt{packmol.inp}, set \texttt{filetype pdb}.
                \item Provide these uniquely named PDB files as input to \texttt{packmol} in the \texttt{structure} blocks.
                \item \texttt{antechamber} must have been run on these same uniquely named PDBs to generate the \texttt{.mol2} files that \texttt{tleap} will use, ensuring atom name consistency.
            \end{itemize}
            \textit{Alternative (if available):} Try a different version or build of \texttt{packmol} that correctly handles \texttt{.mol2} inputs.
        \end{itemize}

    \item \textbf{Issue: \texttt{tleap} Parameter Assignment/Residue Recognition Errors}
        \begin{itemize}
            \item \textit{Symptom:} \texttt{tleap} cannot assign parameters because it doesn't recognize an atom, its type, or the residue it belongs to. Or, it reports splitting a known residue into an unknown part (e.g., "POL and OL").
            \item \textit{Likely Cause(s) \& Solution(s):}
            \begin{itemize}
                \item \textbf{Residue Name Mismatch:} The residue name in \texttt{mixture\_packed.pdb} (from \texttt{packmol}) does not exactly match the name of a unit defined in \texttt{tleap} via \texttt{RES = loadmol2 res.mol2}. Ensure consistency.
                \item \textbf{Atom Name Mismatch:} An atom name in \texttt{mixture\_packed.pdb} (within a specific residue) does not exactly match an atom name within the corresponding \texttt{.mol2} file loaded by \texttt{tleap}. This was why our unique atom naming strategy was important for the Packmol PDB workaround.
                \item \textbf{Problematic Residue Name (Name Conflicts/Splitting):} As encountered with \texttt{POL}, some 3-letter residue names might conflict with internal \texttt{tleap} libraries or parsing rules.
                    \begin{itemize}
                        \item \textit{Solution:} Choose a more unique 3-letter residue name (e.g., we changed \texttt{POL} to \texttt{PRH}). Ensure this new name is used consistently: in the input PDB for \texttt{antechamber}, with the \texttt{-rn} flag in \texttt{antechamber}, as the residue name within the \texttt{.mol2} file, as the residue name in the PDB for \texttt{packmol}, and when defining the unit in \texttt{tleap} (e.g., \texttt{PRH = loadmol2 pol.mol2} where \texttt{pol.mol2} internally defines \texttt{PRH}).
                    \end{itemize}
                \item \textbf{Incorrect \texttt{.mol2} File:} The \texttt{.mol2} file loaded by \texttt{tleap} might be corrupted, not contain proper GAFF2 atom types, or lack charges. Ensure \texttt{antechamber} ran correctly.
                \item \textbf{Missing \texttt{.frcmod} File:} If a \texttt{.frcmod} file was generated but not loaded in \texttt{tleap}, parameters might be missing.
                \item \textbf{\texttt{leaprc.gaff2} not sourced:} If the base force field isn't loaded, no types or parameters will be available.
            \end{itemize}
        \end{itemize}

    \item \textbf{Issue: InterMol/ParmEd GROMACS Topology Issues}
        \begin{itemize}
            \item \textit{Symptom:} \texttt{gmx grompp} fails with "No molecules were defined in the system." Inspection of the \texttt{.top} file shows a single large \texttt{[ moleculetype ]} or missing \texttt{[ molecules ]} section.
            \item \textit{Likely Cause(s):} The conversion tool did not correctly interpret the multiple residue types from the Amber \texttt{prmtop} as separate GROMACS \texttt{moleculetype}s.
            \item \textit{Solution(s) (InterMol CLI worked for us):}
            \begin{itemize}
                \item Use \texttt{python -m intermol.convert --amb\_in mixture.prmtop mixture.inpcrd --gromacs --oname mixture\_intermol}. This invocation correctly generated separate \texttt{[ moleculetype ]} sections for \texttt{MOH} and \texttt{PRH} and the final \texttt{[ system ]} and \texttt{[ molecules ]} sections.
                \item If using ParmEd's Python API or other InterMol API calls, ensure the methods used are appropriate for distinguishing multiple residue types into separate \texttt{moleculetype}s. This can be version-dependent.
            \end{itemize}
        \end{itemize}

    \item \textbf{Issue: GROMACS \texttt{mdrun} Crashes Immediately (LINCS Warnings)}
        \begin{itemize}
            \item \textit{Symptom:} Simulation blows up immediately. LINCS warnings show huge constraint deviations.
            \item \textit{Likely Cause(s):} Severe steric clashes or strained geometry in the starting structure, leading to enormous forces.
            \item \textit{Solution(s):} Perform Energy Minimization. Always run energy minimization (\texttt{gmx grompp -f minim.mdp ...} then \texttt{gmx mdrun -deffnm system\_em}) on the initial structure from \texttt{packmol}/conversion before starting any dynamics (NVT, NPT).
        \end{itemize}

    \item \textbf{Issue: GROMACS NPT Simulation Instability (Barostat Choice)}
        \begin{itemize}
            \item \textit{Symptom:} NPT simulation is unstable, box volume oscillates wildly, then LINCS failures and crash.
            \item \textit{Likely Cause(s):} The system is not yet sufficiently equilibrated (density is far from target), and an "aggressive" barostat like Parrinello-Rahman is used too early.
            \item \textit{Solution(s):} For initial NPT equilibration, use a more stable, "gentler" barostat like Berendsen (\texttt{pcoupl = Berendsen} in \texttt{npt.mdp}). Run with Berendsen until the density stabilizes, then optionally switch to Parrinello-Rahman for longer NPT runs if precise NPT ensemble statistics are needed from that stage.
        \end{itemize}

    \item \textbf{Issue: Gnuplot (or other graphical tools) Display Errors on Headless Server}
        \begin{itemize}
            \item \textit{Symptom:} Plotting commands fail with "could not connect to display" or "cannot open display".
            \item \textit{Likely Cause(s):} The program is trying to open a graphical window, but the SSH session to the headless server does not have X11 display forwarding enabled or working.
            \item \textit{Solution(s):}
            \begin{itemize}
                \item \textbf{Output to File:} Instruct Gnuplot to render to an image file:
                \begin{lstlisting}[language=gnuplot, caption=Gnuplot output to file, numbers=none]
set terminal pngcairo size 800,600
set output 'my_plot.png'
plot 'data.xvg' using 1:2 with lines
unset output # Important to close the file
                \end{lstlisting}
                Then download the \texttt{.png} file using \texttt{scp} or \texttt{sftp}.
                \item \textbf{Enable X11 Forwarding (More complex):} Configure X11 forwarding on the server (\texttt{/etc/ssh/sshd\_config}) and use \texttt{ssh -X} or \texttt{ssh -Y} when connecting from a local Linux/macOS machine (which must have an X server running). For Windows clients, use an X server like Xming or VcXsrv and configure your SSH client (e.g., PuTTY, MobaXterm) for X11 forwarding.
            \end{itemize}
        \end{itemize}

    \item \textbf{Issue: \texttt{scp} / \texttt{sftp} File Transfer Failures}
        \begin{itemize}
            \item \textit{Symptom:} Cannot connect to the server from the local machine to transfer files.
            \item \textit{Likely Cause(s) \& Solution(s):}
            \begin{itemize}
                \item \textbf{Incorrect IP Address or Hostname:} Verify the server's IP/hostname.
                \item \textbf{Server Firewall:} Ensure the server's firewall (e.g., \texttt{ufw} on Ubuntu) allows incoming connections on port 22 (SSH). Command: \texttt{sudo ufw allow 22/tcp}.
                \item \textbf{Local Firewall:} Ensure your local machine's firewall isn't blocking outgoing SSH/SCP connections (less common).
                \item \textbf{SSH Daemon Not Running on Server:} Check \texttt{sudo systemctl status ssh} on the server.
                \item \textbf{Network Issues:} Problems with your local network, router, or internet connection.
                \item \textbf{Using Cloudflare Tunnel (Specific to your case):} If accessing the server via a Cloudflare Tunnel, \texttt{scp} and \texttt{sftp} from your local machine must use the SSH host alias defined in your local SSH config file (e.g., \texttt{scp toricalc:/path/on/server /local/path}), which uses the \texttt{ProxyCommand} to route through \texttt{cloudflared.exe}. Direct connection to the server's private IP will fail if you are not on the same local network.
            \end{itemize}
        \end{itemize}
\end{itemize}

%==================================================================
%   END OF TROUBLESHOOTING
%==================================================================

\section{Conclusion}
\label{sec:conclusion}
% --- START: Placeholder for Conclusion ---
This guide has detailed a comprehensive workflow for simulating the Vapor-Liquid Equilibrium (VLE) of a binary mixture, using methanol and 1-propanol as a practical example. By systematically employing tools from the AmberTools suite (namely \texttt{antechamber}, \texttt{parmchk2}, and \texttt{tleap}) for molecular parameterization with the GAFF2 force field, and \texttt{packmol} for initial system assembly, we prepared a system suitable for molecular dynamics. A crucial conversion step to GROMACS format was robustly handled using InterMol, overcoming common pitfalls associated with multi-component topology generation.

The subsequent GROMACS simulation protocol involved essential equilibration steps: energy minimization to relax initial strains, NVT equilibration to thermalize the system at the target temperature, and NPT equilibration to achieve the correct liquid density under target pressure. The choice of the Berendsen barostat for initial NPT proved vital for stability. Following equilibration, a slab geometry was created using \texttt{gmx editconf}, setting the stage for the NVT VLE production run. Analysis of the VLE trajectory using \texttt{gmx density} allowed for the determination of component density profiles, from which liquid and vapor phase compositions ($x_i, y_i$) can be derived. The equilibrium vapor pressure ($P_{\text{vap}}$) can also be extracted using \texttt{gmx energy}.

The successful execution of this multi-stage process yields a single VLE tie-line ($T, P_{\text{vap}}, x_i, y_i$). By systematically repeating this workflow at various temperatures (or overall compositions), users can generate the necessary data to construct complete VLE phase diagrams (e.g., T-xy, P-xy).

While this workflow is powerful, users should be mindful of the computational expense, especially for long VLE production runs and for generating multiple data points. The accuracy of the results is inherently tied to the quality of the chosen force field (GAFF2). Comparison with experimental data, where available, is always recommended for validation. The troubleshooting section provides guidance for common issues encountered during such complex simulation setups.

Ultimately, this documented methodology provides a robust pathway for computational chemists and chemical engineers to investigate and predict VLE behavior of organic mixtures, contributing valuable data for process design and fundamental thermodynamic understanding.
%==================================================================
%   END OF CONCLUSION
%==================================================================

\appendix
\section{Appendix: Full Scripts and Input Files}
\label{app:scripts}

All scripts and key input files discussed in this guide are available on GitHub. Please refer to the repository for the complete files:
\begin{center}
    \url{https://github.com/Thyco5/gromacs_vle/upload/main} % <<<--- REPLACE THIS WITH YOUR ACTUAL GITHUB LINK
\end{center}
\vspace{1em} % Add some vertical space

For convenience, the content of the main scripts and input files generated throughout this workflow are also listed below.

\subsection*{Example Methanol PDB (\texttt{methanol\_unique.pdb})}
The content for this file is shown in Listing \ref{lst:pdb_moh}.

\subsection*{Example Antechamber Command (Methanol)}
The command used is shown in Listing \ref{lst:antechamber_moh_cmd}.

\subsection*{Example Parmchk2 Command (Methanol)}
The command used is shown in Listing \ref{lst:parmchk2_moh_cmd}.

\subsection*{Example Packmol Input File (\texttt{packmol.inp})}
The content for this file is shown in Listing \ref{lst:packmol_inp}.

\subsection*{Example \texttt{tleap} Input File (\texttt{tleap\_mixture.in})}
The content for this file is shown in Listing \ref{lst:tleap_in}.

\subsection*{Example InterMol Python Script (\texttt{convert\_with\_intermol.py})}
The content for this file is shown in Listing \ref{lst:intermol_py}.
Alternatively, the command line invocation is shown in Listing \ref{lst:run_intermol_cli}.

\subsection*{Example Energy Minimization MDP File (\texttt{minim.mdp})}
The content for this file is shown in Listing \ref{lst:minim_mdp}.

\subsection*{Example NVT Equilibration MDP File (\texttt{nvt\_equil.mdp})}
The content for this file is shown in Listing \ref{lst:nvt_mdp}.

\subsection*{Example NPT Equilibration MDP File (\texttt{npt.mdp})}
The content for this file is shown in Listing \ref{lst:npt_mdp}.

\subsection*{Example VLE Production MDP File (\texttt{vle\_nvt.mdp})}
The content for this file is shown in Listing \ref{lst:vle_nvt_mdp}.

\subsection*{Example Gnuplot Script for Density Profiles}
An example script is shown in Listing \ref{lst:gnuplot_density}.

\subsection*{Example Python Snippet for Average Density Calculation}
A conceptual script is shown in Listing \ref{lst:avg_density_py}.

\subsection*{Example Python Snippet for Average Pressure Calculation}
A conceptual script is shown in Listing \ref{lst:avg_pressure_py}.

% You can add more specific references or full listings here if desired.
% For instance, if you have a master shell script that runs the entire workflow.

\end{document}
```
