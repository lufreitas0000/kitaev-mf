\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{tufte-style-thesis}[]
%\RequirePackage{amsfonts}		% math fonts
%\RequirePackage{amsmath}		% math stuff
%\RequirePackage{mathtools}		% amsmath extension
%\RequirePackage[T1]{fontenc}	% font encoding

%\RequirePackage{fontspec} 
\RequirePackage[utf8]{inputenc}
%\usepackage[utf8]{inputenc}
%\usepackage[latin9]{inputenc}





%%% OPTIONS %%%

\newif\if@justified\@justifiedfalse
\DeclareOption{raggedright}{	% makes all paragraphs left-aligned
	\AtBeginDocument{
		\setlength{\RaggedRightParindent}{1em}
		\RaggedRight
	}
}
\DeclareOption{boxey}{
    \@justifiedtrue
}
\DeclareOption{parskip}{		% changes paragraph indents vor vertical skips
	\AtBeginDocument{\usepackage{parskip}}
	\setlength{\parskip}{\baselineskip}
}
\DeclareOption{casual}{		% removes header numbering --still puts them in toc
	\AtBeginDocument{\setcounter{secnumdepth}{3}}
}
\newif\if@sans\@sansfalse
\DeclareOption{sans}{
	\@sanstrue
	\renewcommand{\familydefault}{\sfdefault}
}
\newif\if@black\@blackfalse
\DeclareOption{black}{
	\@blacktrue
}
\newif\if@colorful\@colorfulfalse
\newif\if@black\@blackfalse
\DeclareOption{colorful}{
	\@colorfultrue
	\@blackfalse
}
\DeclareOption{black}{
	\@blacktrue
	\@colorfulfalse
}

\ExecuteOptions{}
\ProcessOptions\relax
\LoadClassWithOptions{book}

\if@justified
    \AtBeginDocument{
        \justifying
    }
\else
    \AtBeginDocument{
        \setlength{\RaggedRightParindent}{1em}
        \RaggedRight
    }
\fi






% GEOMETRY, TEXTWIDTH
\RequirePackage[a4paper,
    left=20mm,
    top=35mm,
    bottom=35mm,
    headsep=2\baselineskip,
    textwidth=150mm,
    marginparsep=4mm,
    marginparwidth=25mm,
    %textheight=44\baselineskip,
    headheight=2\baselineskip
]{geometry}


%%% PACKAGE DEFINITIONS %%% 

\RequirePackage{amsfonts}		% math fonts
\RequirePackage{amsmath}		% math stuff
\RequirePackage{mathtools}		% amsmath extension
\RequirePackage{emptypage}		% if a page is empty, it is really empty
\RequirePackage{fullwidth}		% for wide environments
\RequirePackage{sidenotes}		% for margin stuff
\RequirePackage[
	hypcap=false				% hypcap=true spits an error
]{caption}						% for caption formatting
\RequirePackage[osf]{libertinus}% main font
\RequirePackage{libertinust1math}
\RequirePackage{gillius}		% sans font
\RequirePackage[
	defaultmono,
	scale=.82
]{droidsansmono}				% mono font
\RequirePackage{ragged2e}		% for better raggedright
\RequirePackage{titlesec}		% header customization
\RequirePackage{titletoc}		% toc customization
\RequirePackage[titles]{tocloft}% lof lot lol customizatios
\RequirePackage{fancyhdr}		% page header customization
\RequirePackage{graphicx}		% for images
\RequirePackage[%				% listings
	% numberbychapter=true,
]{listings}
\RequirePackage[
	protrusion=true,
	expansion=true,
	final,
	tracking,
]{microtype}					% microtypography: fine-tuning in text
%\RequirePackage{physics}		% handy shortcuts for physics
\RequirePackage{xcolor}			% colorz
\RequirePackage[framemethod=TikZ]{mdframed} % for custom code boxes
\RequirePackage{tabularx}		% adaptive columns on tables
%\RequirePackage{booktabs}		% better looking tables
%\RequirePackage{enumitem}		% better looking lists
\RequirePackage{setspace}
\RequirePackage[	hidelinks,
    final,
    breaklinks,
    bookmarks,	pdfusetitle,
]{hyperref}	% automatic references in pdf			
\RequirePackage{etoolbox}		% really idk what this does
\RequirePackage[strict]{changepage}	% needed by sidenotes
\RequirePackage{placeins}		% for floatbarrier
\RequirePackage{xparse}
\RequirePackage{xpatch}

\definecolor{nicepink}{HTML}{ea147e}
\RequirePackage[
	% doi=false,
	% isbn=false,
	% url=false,
%	autocite 	= superscript,
	backend 	= biber,
%	style 		= nature,
%	  citecolor   = nicepink,
	sorting		= none,
	bibencoding = utf8
]{biblatex}
%\RequirePackage{csquotes}
%%% TEST footmisc
%\RequirePackage[    perpage,    symbol*,    flushmargin,    side,    ragged, ]{footmisc}


% \RequirePackage{marginfix}
\makeatletter
  \renewcommand\@makefnmark{%
    \hbox{%
      \ifx\@fnsymbol\FN@fnsymbol@sym
        \expandafter\@firstofone
      \else
        \expandafter\@textsuperscript
      \fi
      {\normalfont\@thefnmark}%
    }%
  }
\makeatother
\RequirePackage{xifthen}


%	C O L O R  D E F I N I T I O N
%
\definecolor{niceRed}{HTML}{980326}
\definecolor{niceBlue}{HTML}{025472}
\definecolor{tooRed}{HTML}{cc0000}
\definecolor{nicepink}{HTML}{ea147e}

\definecolor{darkblue}{cmyk}{100,36,0,40}
\definecolor{darkred}{cmyk}{0, 89, 88, 12}
\definecolor{darkgreen}{cmyk}{56, 0, 58, 31}
\definecolor{darkpink}{cmyk}{0, 91, 46, 8}


\definecolor{pink1}{cmyk}{0,40,31,0}
\definecolor{blue1}{cmyk}{70,20,0,0}
\definecolor{purple1}{cmyk}{19, 71, 0, 5}
\definecolor{orange1}{cmyk}{0, 39 78, 4}
\definecolor{yellow1}{cmyk}{0, 16, 62, 5}
\definecolor{brown1}{cmyk}{0, 35, 71, 50}

\if@colorful
	\renewcommand\thefootnote{\color{nicepink}\arabic{footnote}}
	\renewcommand\thesidenote{\color{nicepink}\arabic{sidenote}}
\else
\fi
\if@black
	\definecolor{niceRed}{HTML}{000000}
	\definecolor{niceBlue}{HTML}{000000}
\else
\fi
\newcommand\red[1]{\textcolor{niceRed}{#1}}






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MARGINS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% formatting sidenote text
\makeatletter
\RenewDocumentCommand\sidenotetext{ o o +m }{%
    \IfNoValueOrEmptyTF{#1}{%
        \@sidenotes@placemarginal{#2}{\setstretch{1.0}\noindent\footnotesize\RaggedRight{\thesidenote{}}~~#3}%
        \refstepcounter{sidenote}%
    }{%
        \@sidenotes@placemarginal{#2}{\noindent\footnotesize\RaggedRight{#1}~~#3}%
    }%
}
\makeatother

\newcommand{\sidetext}[1]
{\marginpar{\setstretch{1.0}\RaggedRight\footnotesize\noindent #1 \vskip 4pt}}

\newenvironment{wide}{
\begin{fullwidth}[outermargin=-\marginparwidth-\marginparsep, width=\linewidth+\marginparwidth+\marginparsep]
}%
{
\end{fullwidth}
}
% sidenotes reset on new
\let\oldchap\chapter
\renewcommand{\chapter}{\setcounter{sidenote}{1}\oldchap}

%%%	FULL ARABIC -- page numbering Arabic numbers  %%%
\makeatletter
\renewcommand{\frontmatter}{\clearpage\@mainmatterfalse}
\renewcommand{\mainmatter}{\clearpage\@mainmattertrue}
%\renewcommand{\frontmatter}{\cleardoublepage\@mainmatterfalse}
%\renewcommand{\mainmatter}{\cleardoublepage\@mainmattertrue}
\makeatother



%%%%%%%%%%%%%%%%%%% PAGE HEADERS %%%%%%%%%%%%%%%%%%%
\renewcommand{\chaptermark}[1]{ \markboth{#1}{} }
\renewcommand{\sectionmark}[1]{ \markright{#1}{} }
\renewcommand{\subsectionmark}[1]{ \markright{#1}{} }
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}
\fancyhead[LO]{	\footnotesize\textls{\MakeUppercase{\sffamily\rightmark}} \hfill \normalfont\normalsize\thepage
}
\fancyhead[RE]{
	\normalfont\normalsize\thepage\hfill \footnotesize\textls{\MakeUppercase{\sffamily\leftmark}}
}
\fancypagestyle{plain}{
\fancyhf{}
\fancyhead[RO]{\thepage}
\fancyhead[LE]{\hfill}
}
\setlength{\headheight}{13.6pt}




%%%%%%%%%%%%%%%%%%%	   TITLE  HEADERS        %%%%%%%%%%%%%%%%%%%

%	default numbering depth
\AtBeginDocument{\setcounter{secnumdepth}{3}}

%	header formatting
\titlespacing*{\part}{.2\linewidth}{.3\textheight}{0pt}
\titlespacing*{\chapter}{0pt}{-20pt}{65pt}
%\titlespacing{\chapter}{0pt}{-60pt}{65pt}

\titleclass{\part}{page}
\titleformat{\part}[display]%[hang]
	{ \huge\raggedleft}
	{\Large\lsstyle\MakeUppercase{\partname}~~\thepart\normalfont}
	{0em}
	{\vskip 12pt \setstretch{0.9}\huge\lsstyle\sffamily\uppercase}


\titleclass{\chapter}{top}
\titleformat{\chapter}
        {\itshape\huge}
        {\normalfont\bfseries\huge\thechapter}
        {1em}
        {\setstretch{0.9}\raggedright\huge\itshape}
\titleformat{\section}
	{\itshape\Large}
	{\normalfont\thesection}
	{1em}
	{\setstretch{0.9}\itshape}
\titleformat{\subsection}
	{\itshape\normalfont}
	{\normalfont\thesubsection}
	{1em}
	{\setstretch{0.9}\itshape}
%\titleformat{\subsubsection}[runin]
%	{}
%	{\itshape\thesubsubsection}
%	{1ex}
%	{\itshape}[~~$\cdot$]

\titleclass{\subsubsection}{straight}
\titleformat{\subsubsection}%
  [hang]% shape
  {\normalfont\large\itshape}% format applied to label+text
  {\thesubsubsection}% label
  {1em}% horizontal separation between label and title body
  {}% before the title body
  []% after the title body



%%%       	MAKE TITLE              %%%

% 	custom info commands like @title, @author...
% 	info commands definition
\makeatletter 
% 	subtitle --title already exists
\newcommand\subtitle[1]{\renewcommand\@subtitle{#1}}
\newcommand\@subtitle{}
% 	university
\newcommand\university[1]{\renewcommand\@university{#1}}
\newcommand\@university{}
% 	lab
\newcommand\lab[1]{\renewcommand\@lab{#1}}
\newcommand\@lab{}

% 	supervisors, jury members, etc.
\newcommand\people{}
\def\zfoo#1\relax{\long\def\people##1{#1}}
\newcommand\updatepeople[1]{%
   \expandafter\zfoo\people{##1}##1#1\relax
}
\newcommand\person[3]{
    \updatepeople{#1 & #2 & #3 \\}
}
%%%%%%%%%%%%%%%%%%%%%%  Insert Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	logos  
\newcommand\logos{}
\def\xfoo#1\relax{\long\def\logos##1{#1}}
\newcommand\updatelogos[1]{%
   \expandafter\xfoo\logos{##1}##1#1\relax
}
\newcommand\logo[1]{
    \updatelogos{%\hspace{20pt}\includegraphics[height=50pt]{#1}
    \hspace*{-20pt}\includegraphics[width=0.45\linewidth]{#1}\hfill
    }
}

% 	Doctoral thesis (subsubtitle)
\newcommand\type[1]{\renewcommand\@type{#1}}
\newcommand\@type{Doctoral thesis}

% 	shoutouts
\newcommand\shoutouts[1]{\renewcommand\@shoutouts{#1}}
\newcommand\@shoutouts{}



%%%%%%%%%%%%%%%% 	ACTUAL MAKE TITLE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand\@maketitlebis{%
\makeatletter
\setstretch{1.0}
\begin{titlepage}
%\setcounter{page}{1}
    \hfill
    \begin{minipage}{0.99\linewidth}
	\raggedleft
	{\large\scshape
		\MakeLowercase{\@university}\\[5pt]
		\MakeLowercase{\@lab}%
        }%
        
        {    \vspace{5mm}
            \logos{}%  inset uni logo
            \vspace*{-40mm}
        }
        
	\vskip 80pt
	{\Large\itshape\@author}
	\vskip 40 pt
	{%
			\huge
			\if@colorful
				\color{niceRed}
			\else
			\fi
			\@title\par
			}%
	\vskip 8pt
	{\Large\textbf\@subtitle\par}%
	\vskip 20pt
	{\Large \@type}%
	\vskip 40pt
	{\itshape\@date}
		\end{minipage}
		\vfill
		\raggedleft
		\begin{tabular}{rll}
			\people{}
		\end{tabular}
		\vskip 40pt
		%\logos{}
	\end{titlepage}
	\if@colorful
		\hypersetup{
			%colorlinks = true,
			linkcolor = nicepink,
			anchorcolor = nicepink,
                linktoc=all,  
			urlcolor = nicepink,
                pdftitle={de Freitas' Thesis},
                pdfauthor={L. de Freitas}
                pdfpagemode=FullScreen, 
                filecolor=cyan,      
                urlcolor=RoyalBlue,
                citecolor=Purple,
		}
	\fi
 
	% copyright
	\thispagestyle{empty}
	\hfill
	\vfill
	\makeatletter
	\noindent\begin{flushleft}\@author,\\ \textit{\@title,}\\ \@subtitle,\\ \@date.
    \end{flushleft}
	\makeatother

	% dedication
	%\cleardoublepage
	%\thispagestyle{empty}
	%~\vskip .3\textheight
	%\begin{flushright}
	%\large\itshape\@shoutouts
	%\end{flushright}

	\setstretch{1.12}
%	\makeatother
}
\renewcommand\maketitle{
	\@maketitlebis
	\thispagestyle{empty}
}
%\makeatother


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	COLOPHON
\AtEndDocument{
%	\cleardoublepage
	%\pagestyle{empty}
	%\footnotesize\raggedleft
	%\hfill
	%\vfill
	%\noindent EDIT COLOPHON %This document was typeset using \LaTeX{} and the \texttt{tufte-style-thesis} class. 
	%\vskip 8pt
	%\centering \url{ }.\\
	%\raggedleft\noindent   
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	F I G U R E S  A N D  S T U F F
%
% 	caption formatting
\DeclareCaptionLabelFormat{lc}{\MakeLowercase{#1}~#2}
\DeclareCaptionLabelSeparator{mysep}{.\quad}
\DeclareCaptionLabelSeparator{mysepside}{.\\}
\captionsetup{
	labelfont = {bf},
	%labelformat = lc,
	labelsep = mysep,
	justification=RaggedRight,
	font={footnotesize, stretch=1},
	singlelinecheck=on,
}
\DeclareCaptionStyle{sidecaption}{
    labelfont = {bf},
	%labelformat = lc,
	labelsep = mysepside,
	justification=RaggedRight,
	font={footnotesize, stretch=1},
	singlelinecheck=off,
}
\DeclareCaptionStyle{marginfigure}{
    labelfont = {bf},
	%labelformat = lc,
	labelsep = mysepside,
	justification=RaggedRight,
	font={footnotesize, stretch=1},
	singlelinecheck=off,
}
\DeclareCaptionStyle{widefigure}{
    labelfont = {bf},
	%labelformat = lc,
	labelsep = mysep,
	justification=RaggedRight,
	font={footnotesize, stretch=1},
	singlelinecheck=off,
}
\DeclareCaptionStyle{margintable}{
    labelfont = {bf},
	%labelformat = lc,
	labelsep = mysep,
	justification=RaggedRight,
	font={footnotesize, stretch=1},
	singlelinecheck=off,
}
\DeclareCaptionStyle{widetable}{
    labelfont = {bf},
	%labelformat = lc,
	labelsep = mysep,
	justification=RaggedRight,
	font={footnotesize, stretch=1},
	singlelinecheck=off,
}
%	figure macros
\AtBeginDocument{\usepackage[export]{adjustbox}} % for 'outer' mode
\newcommand{\textfig}[4][]{%
	\FloatBarrier
	\begin{figure}[!htb]%
		\sidecaption{#3\label{#4}}%
		\includegraphics[width = #1\linewidth, outer]{#2}%
	\end{figure}%
}
\newcommand{\marginfig}[3]{%
	\FloatBarrier%
	\begin{marginfigure}%
		\includegraphics[width = \linewidth]{#1}%
		\caption{#2\label{#3}}%
	\end{marginfigure}%
}
\newcommand{\widefig}[4][]{%
	\FloatBarrier%
	\begin{figure*}[!htb]%
		\includegraphics[width=#1\linewidth, outer]{#2}%
		\sidecaption{#3\label{#4}}%
	\end{figure*}%
}
\newcommand{\plainfig}[4][]{%
	\FloatBarrier%
	\begin{figure}[!htb]%
		\includegraphics[width=#1\linewidth, left]{#2}%
		\caption{#3\label{#4}}%
	\end{figure}%
}
\newcommand{\plainwidefig}[4][]{%
	\FloatBarrier%
	\begin{figure*}[!htb]%
		\includegraphics[width=#1\linewidth, left]{#2}%
		\caption{#3\label{#4}}%
	\end{figure*}%
}
%	table macros
\newcolumntype{Y}{>{\centering\arraybackslash}X} % tabularx column type
\newenvironment{texttable}[2]{%
	\begingroup%
	\small%
	\marginpar{%
		\vskip 8pt%
		\captionof{table}{#1}%
		\label{#2}%
	}

	\noindent%
}
{
	\endgroup%
}


%%%%%%%%%%%%%%%%%%%%%%             TOC             %%%%%%%%%%%%%%%%%%%%%%
\let\oldtoc\tableofcontents
\renewcommand{\tableofcontents}{%
	\begin{quote}\setstretch{0.5}
	\oldtoc
	\end{quote}
}
%%%%%%%%%%%%%%%%%%%%%             TOC depth             %%%%%%%%%%%%%%%%%%%%%
\AtBeginDocument{\setcounter{tocdepth}{0}} %0=chapters, 1=sections, ...
\renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}%


%	TOC  C U S T O M I Z A T I O N
\contentsmargin[1cm]{0cm}
\titlecontents*{part}[0em]
	{\large\scshape}
	{\scshape\thecontentslabel}
	{\vskip 8pt\large\scshape\lowercase}
	{}[\vskip 4pt]
\titlecontents{chapter}[0em]
	{\large\vspace{6pt}}
	{\makebox[-1em][r]{\thecontentslabel}\large\hspace{1em}}
	{\large}
	{~~~~\textcolor{nicepink}{\contentspage}}[\vskip 1pt]
\titlecontents*{section}[0em]
	{\normalsize}%it was small before
        {\thecontentslabel~~}
	{\small}
	{~~\dotfill~~\textcolor{nicepink}{\contentspage}}
	[\vskip 1pt]
\titlecontents{subsection}[0em]
 	{\vspace{4pt}}
 	{\normalfont{\thecontentslabel}\enspace\itshape}
 	{\itshape}
 	{\qquad\contentspage}[\vskip 0pt]
% \titlecontents*{subsection}[0em]
% 	{\vskip 0pt\small}
% 	{{\thecontentslabel}~~\itshape}
% 	{\itshape}
% 	{~~{\contentspage}}
% 	[\qquad]
 \titlecontents*{subsubsection}[1em]
 	{\vspace{0pt}}
 	{\normalfont\thecontentslabel\enspace\itshape}
 	{\itshape}
 	{\quad\quad\contentspage}[\vskip 0pt]


% TOCLOFT (classicthesis)
\renewcommand{\cftfigpresnum}{\scshape\MakeTextLowercase}%
\renewcommand{\cftfigfont}{\normalfont}%
\renewcommand{\cftfigleader}{\hspace{2.5em}}%
\renewcommand{\cftfigafterpnum}{\cftparfillskip}%
\renewcommand{\cftfigpresnum}{\figurename~}%Fig.~}
\newlength{\figurelabelwidth}
\settowidth{\figurelabelwidth}{\cftfigpresnum~999}
\addtolength{\figurelabelwidth}{1em}
\cftsetindents{figure}{0em}{\figurelabelwidth}
% tables
\renewcommand{\cfttabfont}{\normalfont}%
\renewcommand{\cfttableader}{\hspace{2.5em}}%
\renewcommand{\cfttabafterpnum}{\cftparfillskip}%
\renewcommand{\cfttabpresnum}{\tablename~}%Tab.~}
\newlength{\tablelabelwidth}
\settowidth{\tablelabelwidth}{\cfttabpresnum~999}
\addtolength{\tablelabelwidth}{1em}
%\cftsetindents{table}{0em}{\tablelabelwidth}
\cftsetindents{table}{0em}{\figurelabelwidth}
% listings
\newcommand\listingsname{Listing}
\renewcommand\lstlistingname{Listing}
\renewcommand\lstlistlistingname{List of Listings}
\newlistof{listings}{lol}{\lstlistlistingname}%
\renewcommand{\cftlistingspresnum}{\scshape\MakeTextLowercase}%
\renewcommand{\cftlistingsfont}{\normalfont}%
\renewcommand{\cftlistingspresnum}{\lstlistingname~}%
\renewcommand{\cftlistingspagefont}{\normalfont}%
\renewcommand{\cftlistingsleader}{\hspace{1.5em}}%
\renewcommand{\cftlistingsafterpnum}{\cftparfillskip}%
\newlength{\listingslabelwidth}%
\settowidth{\listingslabelwidth}{\cftlistingspresnum~999}%
\addtolength{\listingslabelwidth}{2.5em}%
\cftsetindents{listings}{0em}{\figurelabelwidth}%
% \let\l@lstlisting\l@listings%
\let\lstlistoflistings\listoflistings%
% blue page numbers
\renewcommand{\cftfigpagefont}{\color{nicepink}}
\renewcommand{\cfttabpagefont}{\color{nicepink}}
\renewcommand{\cftlistingspagefont}{\color{nicepink}}
% TOCLOFT over










%	B I B L I O G R A P H I C
%	      R E F E R E N C E S 
\makeatletter
\newcommand\sidecite[2][]{\autocite{#2}\sidetext{\ifthenelse{\equal{#1}{}}{}{\vskip #1}%
\@for\cle:={#2}\do{\makebox[0pt][r]{\autocite{\cle}}\fullcite{\cle}.\\\vskip 4pt}%
}}
\newcommand\sidetextcite[2][]{\sidetext{\ifthenelse{\equal{#1}{}}{}{\vskip #1}%
\@for\cle:={#2}\do{\makebox[0pt][r]{\autocite{\cle}}\fullcite{\cle}.\\\vskip 4pt}%
}}
\makeatother

\DeclareCiteCommand{\fullcite}
  {\usebibmacro{prenote}}
  {\clearfield{url}%
   \clearfield{pagetotal}%
   \clearfield{edition}%
   \usedriver
     {\DeclareNameAlias{sortname}{default}}
     {\thefield{entrytype}}}
  {\multicitedelim}
  {\usebibmacro{postnote}}


%	flush       T I G H T  L I S T S
%
\AtBeginDocument{
	\setlist{topsep=0pt, itemsep=0pt, noitemsep}
}
% 	diamond bullets
\renewcommand{\labelitemi}{$\vcenter{\hbox{\small$\diamond$}}$}

%	CODE
%
%   colors
\definecolor{bg}{rgb}{.95,.95,.95}
\definecolor{key}{HTML}{333399}
\definecolor{comment}{HTML}{6b6b99}
\definecolor{codegray}{rgb}{0.4,0.4,0.4}
\definecolor{string}{HTML}{2b4f2b}
%   boxes
\mdfdefinestyle{code_box}{%
    roundcorner= 2pt,
	backgroundcolor=white,
    linecolor = lightgray,
    linewidth = .5pt,
    topline = true,
    bottomline = true,
    leftline = true,
    rightline = true,
	innertopmargin=.2em,
	innerbottommargin=.2em,
}
%   listing styles
\lstdefinestyle{unnumbered}{
    commentstyle=\color{comment},
    keywordstyle=\color{key},
    stringstyle=\color{string},
    basicstyle=\ttfamily\normalsize\lst@ifdisplaystyle\small\fi,
    breakatwhitespace=false,
    breaklines=true,
    keepspaces=true,
    showspaces=false,
    showstringspaces=false,
    showtabs=false,
    tabsize=4,
    framexleftmargin=30pt,
}
\lstdefinestyle{numbered}{
    commentstyle=\color{comment},
    keywordstyle=\color{key},
    stringstyle=\color{string},
    basicstyle=\ttfamily\normalsize\lst@ifdisplaystyle\small\fi,
    breakatwhitespace=false,
    breaklines=true,
    keepspaces=true,
    numbers=left,
    numberstyle = \color{codegray}\footnotesize\ttfamily,
    numbersep=10pt,
    showspaces=false,
    showstringspaces=false,
    showtabs=false,
    tabsize=4,
    framexleftmargin=30pt,
    framextopmargin=100pt,
    framexbottommargin=30pt,
}
%   code macros
\lstnewenvironment{codebox}[1]{%
	\FloatBarrier
	\mdframed[style = code_box, innerleftmargin=10pt]%
	\lstset{style=unnumbered, language=#1}
}
{
    \endmdframed
}
\lstnewenvironment{altcodebox}[1]{%
	\FloatBarrier
	\mdframed[style = code_box, innerleftmargin=10pt]%
	\lstset{style=unnumbered, language=#1}
}
{
    \endmdframed
}
\lstnewenvironment{codeboxnum}[1]{%
	\FloatBarrier
	\mdframed[style = code_box, innerleftmargin=25pt]%
	\lstset{style=numbered, language=#1}
}
{
    \endmdframed
}
\newcommand{\inlinecode}[1]{\lstinline[style=unnumbered, language=#1]}
%	snippet environments
\lstnewenvironment{snippet}[3]{%
	\marginpar{~\vskip 10pt\captionof{listings}{#2\label{#3}}}
	\mdframed[style = code_box, innerleftmargin=10pt]%
	\lstset{style=unnumbered, language=#1}
}
{
    \endmdframed
}
\lstnewenvironment{snippetnum}[3]{%
	\marginpar{~\vskip 10pt\captionof{listings}{#2\label{#3}}}
	\mdframed[style = code_box, innerleftmargin=25pt]%
	\lstset{style=numbered, language=#1}
}
{
    \endmdframed
}


%	I N - T E X T  F I N E T U N I N G
%
% widows and orphans
\clubpenalty = 10000
\widowpenalty = 10000
\displaywidowpenalty = 10000


%	M I S C E L L A N E O U S
%
\providecommand{\toight}{%
	\tightlist%
}
%	minipages with margin notes
\usepackage{minipage-marginpar}
%	no page break environment
\newenvironment{absolutelynopagebreak}
  {\par\nobreak\vfil\penalty0\vfilneg
   \vtop\bgroup}
  {\par\xdef\tpd{\the\prevdepth}\egroup
   \prevdepth=\tpd}

%	M A T H  &  P H Y S I C S
%	S H O R T C U T S
%\iffalse
% uncomment this \iffalse
%and the \fi at the end to disable (or delete)

%\renewcommand{\phi}{\varphi}
%\renewcommand{\epsilon}{\varepsilon}
\newcommand{\N}{\mathbb{N}}     % shortcuts for usual number sets
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\C}{\mathbb{C}}
% powers of 10%\newcommand{\ee}[1]{\cdot 10^{#1}} 
\newcommand{\vect}[1]{\:\overrightarrow{#1}}
%fi

% TODO
\newlistof{todos}{tod}{\bfseries \textls{TO DO} list}
\let\oldlistoftodos\listoftodos
\renewcommand{\listoftodos}{\bgroup\color{tooRed}{\oldlistoftodos}\egroup}
\newcommand{\todolist}{\listoftodos}

\newcounter{todo}
\setcounter{todo}{1}
\newcommand{\todo}[1]{
    \addcontentsline{tod}{todos}{\protect\numberline{\bfseries\#\thetodo}#1}
    \bgroup\color{tooRed}
    \noindent\large\textsc{todo \#\arabic{todo}}~~\hrulefill%
    \refstepcounter{todo}%
    \smallskip

    \normalsize\noindent #1

    \noindent\hrulefill
    \egroup
}



%		STRETCH
%
\AtBeginDocument{%
	\setstretch{1.15}
}


