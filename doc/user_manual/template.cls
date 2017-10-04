\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mathese}[]
% version 1.0 Olivier Pothier 1999/09/15
% version 1.1 Catherine Lamy 2000/04/18
% translation to english : Charles Rocabert 2016/09/25

%---------------------------------------
%---- on se base sur la classe book ----
%---- sur papier a4, 2 faces, 1 col ----  
%---------------------------------------
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions
\LoadClass{book}

%----------------------
%---- les packages ----
%----------------------
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[francais,english]{babel}

\usepackage{tocloft}
\setlength{\cftfignumwidth}{4em}
%\usepackage[rightcaption]{sidecap}
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}

\usepackage[english]{algorithm2e}
\usepackage{epigraph}
\usepackage{array}
\usepackage{placeins}
\usepackage{graphics,graphicx,epsfig}
\usepackage{enumerate}
\usepackage[usenames,dvipsnames]{color}
\usepackage{natbib}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{amsthm}
\usepackage{euscript}
\usepackage{battail} 
\usepackage{fancyhdr}
\usepackage{latexsym}
\usepackage{array,tabularx}
\usepackage{float}
\usepackage{url}
\usepackage{footnpag}
\usepackage{wrapfig}
\usepackage{multirow}
\usepackage[margin=0pt,font={small,sf},labelfont=bf]{caption}
%\usepackage[margin=0pt,font={small,sf},labelfont=bf]{caption}
\usepackage[nottoc, notlof, notlot]{tocbibind}
\usepackage{setspace} 
\usepackage{pdflscape}
\usepackage{subfig} 
\usepackage{geometry}
\usepackage{changepage}
\usepackage{lipsum}
\usepackage{array}
\usepackage{nicefrac}
\usepackage{bm}
\usepackage{fancybox}
%\usepackage{footbib}
\usepackage[table]{xcolor}
\newcolumntype{C}[1]{>{\centering\arraybackslash}m{#1}}
\usepackage{eurosym}

\usepackage{tocstyle}
\usetocstyle{standard}

%%%%%%%%%%%%%%%%%%%%%%%%% nomenclature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usepackage[english,refpage]{nomencl}
\renewcommand{\pagedeclaration}[1]{\unskip\dotfill{#1}}
\usepackage{lipsum}
\newcounter{symbol_counter}
\newcommand{\nomentry}[2]{% to sort by pages
  \nomenclature[\two@digits{\value{symbol_counter}}]{#1}{#2}%
  \addtocounter{symbol_counter}{1}
} % Ref: http://tex.stackexchange.com/questions/30930/how-to-output-a-counter-with-leading-zeros
%\newcommand{\nomentry}[2]{% to sort normally
%  \nomenclature{#1}{#2}%
%}
\renewcommand{\nompreamble}{This page contains the list of all notations used in this manuscript, with the page number of their first appearance.}
\makenomenclature

%%%%%%%%%%%%%%%%%%%%%%% user defined commands %%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\fracinline}[2]{\raisebox{0.4ex}{$#1$} / \raisebox{-0.7ex}{$#2$}}

\theoremstyle{plain}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}

\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{property}[theorem]{Property}

\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------
%---- la taille de la page ----
%------------------------------
%\voffset 0.0cm
%\hoffset 0.0cm
\textheight 23.5cm
\textwidth 16cm
\topmargin 0.0cm
\oddsidemargin 0cm
\evensidemargin 0cm


% Alter some LaTeX defaults for better treatment of figures:
    % See p.105 of "TeX Unbound" for suggested values.
    % See pp. 199-200 of Lamport's "LaTeX" book for details.
    %   General parameters, for ALL pages:
    \renewcommand{\topfraction}{0.9}	% max fraction of floats at top
    \renewcommand{\bottomfraction}{0.8}	% max fraction of floats at bottom
    %   Parameters for TEXT pages (not float pages):
    \setcounter{topnumber}{2}
    \setcounter{bottomnumber}{2}
    \setcounter{totalnumber}{4}     % 2 may work better
    \setcounter{dbltopnumber}{2}    % for 2-column pages
    \renewcommand{\dbltopfraction}{0.9}	% fit big float above 2-col. text
    \renewcommand{\textfraction}{0.07}	% allow minimal text w. figs
    %   Parameters for FLOAT pages (not text pages):
    \renewcommand{\floatpagefraction}{0.7}	% require fuller float pages
	% N.B.: floatpagefraction MUST be less than topfraction !!
    \renewcommand{\dblfloatpagefraction}{0.7}	% require fuller float pages




%----------------------------------------------------------------
%---- pour la sombre histoire que \chapter* ne fait pas appel
%---- a chaptermark, on a toujours \leftmark qui est en 
%---- majuscules (chiant au possible, on remarque au passage
%---- que T. Oetiker, dans ``The Not So Short Introduction to
%---- Latex 2e'' a laisse tomber...). On definit donc un ersazt
%---- de \chapter*, \chapternonum, qui fait appel a \markboth
%---- avec deux arguments identiques (notons au passage que 
%---- \chapter, lui, comme un con, ne passe qu'un argument, le
%---- premier, a markboth au travers de \chaptermark. Il faudrait
%---- aussi reecrire \chapter, demain peut-etre :-))
%----------------------------------------------------------------
\newcommand{\chapternonum}[1]{\chapter*{#1\markboth{#1}{#1}}}

\renewcommand\thepart{\Alph{part}}
\renewcommand{\thechapter}{\Roman{chapter}} 
\renewcommand{\thesection}{\thechapter.\arabic{section}} 
\renewcommand{\thesubsection}{\thesection.\arabic{subsection}} 

%-----------------------------------
%---- on se charge des headings ----
%-----------------------------------
\pagestyle{fancy}
%\bibliographystyle{fancy}
%\bibliographystyle{apalike}
%---- on doit redefinir \chaptermark et \sectionmark pour qu'ils
%---- ne mettent pas tout en majuscule. On vire au passage tous les
%---- tests debiles fait dans leur version originale, et on corrige
%---- le fait que \chaptermark ne passait qu'un seul argument a
%---- markboth, ce qui posait probleme dans certains cas, rares il est
%---- vrai.
%---- Je ne veux pas afficher \chaptername dans le header, et on
%---- laisse tout dans la casse voulue par l'auteur... L'argument
%---- passe par \chapter a \chaptermark est son argument propre (a
%---- savoir le nom du chapitre. On fait de meme pour section.
\renewcommand{\chaptermark}[1]{\markboth {\thechapter.\ #1}{\thechapter.\ #1}}
\renewcommand{\sectionmark}[1]{\markright{\thesection. \ #1}}
\fancyhf{}
%----       le num"ro de la page       ----
%---- a gauche sur les pages de gauche ----
%---- a droite sur les pages de droite ----
\fancyhead[LE,RO]{\thepage}
%---- \leftmark, cad le chapitre, a droite des pages de gauche
\fancyhead[RE]{\textsc{\leftmark}}
%---- \rightmark, cad la section, a gauche des pages de droite
\fancyhead[LO]{\textsc{\rightmark}}
%---- les initiales ----%
%\fancyfoot[RE,LO]{\footnotesize C.K.}
%---- la version ----%
%\fancyfoot[CE,CO]{\footnotesize version 1.0}
%---- la date ----%
%\fancyfoot[LE,RO]{\footnotesize \today}
%---- les barres ----%
\renewcommand{\headrulewidth}{1pt}
%\renewcommand{\footrulewidth}{1pt}
\addtolength{\headheight}{2.5pt}
%---- on redefinit le style plain (\chapter le force pour la
%---- premiere page de tout chapitre
\fancypagestyle{plain}{
  \fancyhead[RE]{}
  \fancyhead[LO]{}
}


\setlength{\parindent}{0pt}
\setlength{\parskip}{1.2em}

%--------------------------------------------
%---- definitions de nouvelles commandes ----
%--------------------------------------------
\let\urlorig\url
\renewcommand{\url}[1]{%
   \begin{otherlanguage}{english}\urlorig{#1}\end{otherlanguage}%
}

\DeclareGraphicsExtensions{.eps,.png,.pdf}

% Pour que la page paire, ajoutee afin de commencer les chapitres
% sur une page impaire, soit sans header/footer
\newcommand{\clearemptydoublepage}{\newpage{\pagestyle{empty}%
                                   \cleardoublepage}}

\newenvironment{myList}%
{ \begin{list}%
	{$\bullet$}%
	{\setlength{\labelwidth}{30pt}%
	 \setlength{\leftmargin}{30pt}%
	 \setlength{\itemsep}{\parsep}}}%    
{ \end{list} }


\newenvironment{myQuote}[2]%
    {\begin{list}{}{\leftmargin#1\rightmargin#2}\item{}}%
    {\end{list}}


% pour faire des lignes de tableau + epaisses
%%%% debut macro %%%%
\makeatletter
\def\hlinewd#1{%
\noalign{\ifnum0=`}\fi\hrule \@height #1 %
\futurelet\reserved@a\@xhline}
\makeatother
%%%% fin macro %%%%

%% Pour changer les marges localement
\newenvironment{changemargin}[2]{\begin{list}{}{%
\setlength{\topsep}{0pt}%
\setlength{\leftmargin}{0pt}%
\setlength{\rightmargin}{0pt}%
\setlength{\topmargin}{0pt}%
\setlength{\listparindent}{\parindent}%
\setlength{\itemindent}{\parindent}%
\setlength{\parsep}{0pt plus 1pt}%
\addtolength{\leftmargin}{#1}%
\addtolength{\rightmargin}{#2}%
%\addtolength{\topsep}{#2}%
}\item }{\end{list}}
