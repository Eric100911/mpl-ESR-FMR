 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % thuemp.cls
  % thuemp: Experiments of Modern Physics in Tsinghua University

  % version 1.1  (2022年12月)

  % 清华大学近代物理实验实验报告LaTeX模板（非官方）
  % LaTeX Template (unofficial) for Experiments in Modern Physics of Tsinghua University

  % by Mingyu Li
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Version 1.1 [2022.12.30]: add English header part, including title, authors, address, abstract and keyword in English

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1) 参考了《控制理论与应用》提供的LATEX模板  http://jcta.alljournals.ac.cn/uploadfile/cta_cn/20170419/kzllyy%20template20170419-2.9.zip
  % 2) BHOSC   BUAAthesis  https://github.com/BHOSC/BUAAthesis/
  % 3) 北航学报 http://bhxb.buaa.edu.cn/UserFiles/File/%E5%8C%97%E8%88%AA%E5%AD%A6%E6%8A%A5%E6%A8%A1%E6%9D%BF17.1.16(1).doc
  % 4) 北航学报自然科学版LaTeX模板（非官方） https://github.com/Htallone/JBUAA
  % 5) 老师提供的 清华近代物理实验 实验报告模版word版本(2021春版)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{thuemp}
              [2021/03/28 v1.0
The LaTeX template for Experiment Report in Modern Physics of Tsinghua University]
\typeout{Document Class `thuemp' v1.0 by Mingyu Li (2021/03)}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{ctexart}}
\ProcessOptions\relax

\LoadClass[twocolumn, twoside, a4paper, UTF8]{ctexart}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usepackage{fancyhdr}
\usepackage{amsmath,amsfonts,amssymb,graphicx}    % EPS 图片支持
\usepackage{unicode-math}
\usepackage{subfigure}   % 使用子图形
\usepackage{indentfirst} % 中文段落首行缩进
\usepackage{bm}          % 公式中的粗体字符(用命令\boldsymbol）
\usepackage{indentfirst} % 中文首段缩进
\usepackage{abstract}    % 2栏文档，一栏摘要及关键字宏包
\usepackage{amsthm}      % 使用定理
\usepackage{booktabs}    % 使用表格
\usepackage{siunitx}
\usepackage{tikz}
\usepackage{titlesec}
\usepackage{times}
\usepackage{wasysym}
\usepackage{pifont}
\usepackage{ccaption}
\usepackage{float}
\usepackage{calc}
\usetikzlibrary{calc, shapes,arrows,positioning}
\usetikzlibrary{circuits.ee.IEC}
% 在导言区添加 circuitikz 包
\usepackage{circuitikz}
\usepackage{environ}
\usepackage{lmodern}
\usepackage{unicode-math}
\usepackage[sort&compress]{gbt7714}
\usepackage{anyfontsize}
\usepackage{hyperref}
\usepackage{tabu}
\usepackage{tabularx}
\usepackage{multirow}
\usepackage{multicol}
\usepackage{longtable}
\usepackage{makecell}


\setmainfont{Times New Roman}

\newcommand{\erhao}{\fontsize{22pt}{\baselineskip}\selectfont}      % 字号设置
\newcommand{\sihao}{\fontsize{14pt}{\baselineskip}\selectfont}      % 字号设置
\newcommand{\xiaosihao}{\fontsize{12pt}{\baselineskip}\selectfont}  % 字号设置
\newcommand{\wuhao}{\fontsize{10.5pt}{\baselineskip}\selectfont}    % 字号设置
\newcommand{\xiaowuhao}{\fontsize{9pt}{\baselineskip}\selectfont}   % 字号设置
\newcommand{\liuhao}{\fontsize{7.875pt}{\baselineskip}\selectfont}  % 字号设置
\newcommand{\myxiaoliuhao}{\fontsize{6.5pt}{\baselineskip}\selectfont}  % 字号设置
\newcommand{\qihao}{\fontsize{5.25pt}{\baselineskip}\selectfont}    % 字号设置

\renewcommand{\tablename}{表}
\renewcommand{\figurename}{图}
\captiondelim{~~}
\setlength{\abovecaptionskip}{2pt}
\setlength{\belowcaptionskip}{2pt}

\DeclareRobustCommand{\ttfamily}{\fontencoding{T1}\fontfamily{lmtt}\selectfont}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    重定义页面边距，使其符合中文刊物习惯。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usepackage{lettrine}
\setlength{\parindent}{2em}
\setlength{\voffset}{-4mm}
\setlength{\parskip}{1pt}
\renewcommand{\baselinestretch}{1.15} %定义行间距1.15
\setcounter{page}{1} \makeatletter
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {1.5ex }%
                                   {1.5ex}%
                                   {\normalfont\bfseries\sihao}}
\renewcommand\subsection{\@startsection {subsection}{1}{\z@}%
                                   {-1pt }%
                                   {1pt}%
                                   {\normalfont\bfseries\heiti\wuhao}}
\renewcommand\subsubsection{\@startsection {subsubsection}{1}{\z@}%
                                   {-1pt }%
                                   {1pt}%
                                   {\normalfont\wuhao}}



\usepackage{geometry}
\newgeometry{
  top=25.4mm, bottom=22.5mm, left=20mm, right=20mm,
  headsep=5mm, headheight=10mm, footskip=10mm,
}
\savegeometry{mastergeometry}
\loadgeometry{mastergeometry}

\setlength{\columnsep}{10mm}


\newcommand{\positiontextbox}[4][]{
  \begin{tikzpicture}[remember picture,overlay]
    \node[inner sep=0pt,right,align=left, #1] at ($(current page.north west) + (#2,-#3)$) {#4};
  \end{tikzpicture}
}

\ctikzset{bipoles/resistor/height=0.3} % 调整电阻符号大小
\ctikzset{bipoles/resistor/width=0.5}

\NewEnviron{empAbstract}{
\vspace{-3em}
\begin{center}
\parbox{\textwidth}{
\wuhao\setlength{\parindent}{2em}
\newlength{\somelength}

{\heiti\wuhao \makebox[\widthof{【关键词】}][s]{【摘 要】}}
{\xiaowuhao 
\BODY
}

{\heiti\wuhao \makebox[\widthof{【关键词】}][s]{【关键词】}} \heiti \wuhao \empKeyword

\vspace*{1.3cm}}
\end{center}
}

\NewEnviron{empAbstractEn}{
\begin{center}
\parbox{\textwidth}{
\wuhao\setlength{\parindent}{2em}

\textbf{\wuhao Abstract:}
{\xiaowuhao 
\BODY
}

\textbf{\wuhao Key Words:} \textbf{ \wuhao \empKeywordEn}

\vspace*{1.3cm}}
\end{center}
}

\def\empKeyword{}
\newcommand{\Keyword}[1]{
    \def\empKeyword{#1}
}

\def\empKeywordEn{}
\newcommand{\KeywordEn}[1]{
    \def\empKeywordEn{#1}
}


%%题目
\newcommand\emptitle[1]{
\title{
\vspace{0.5cm} \erhao #1 \vspace{-0.2cm}
}
}

%%作者
\newcommand\empauthor[2]{
\author{
\sihao\kaishu #1\makebox{$^{\text{*}}$} ， #2\makebox{{$^\dagger$}} \\[0.1cm]
\liuhao \kaishu ~~清华大学~~致理书院，北京~~100084\ 
}
}

%%英文题目
\newcommand\emptitleEn[1]{
\begin{center}
{
\vspace{0.5cm} \sihao #1 \vspace{-0.2cm}
}
\end{center}
}

%%英文作者
\newcommand\empauthorEn[2]{
\begin{center}
{
\centering
\textbf{\wuhao #1 ， #2} \\
\liuhao Department of Physics, Tsinghua University, Beijing 100084, China\ 
}
\end{center}
}

\date{}

%%首页脚注
\newcommand\empfirstfoot[4]{
\positiontextbox{2.0cm}{27cm}{
\noindent\rule{4cm}{.5pt}\\[0.5ex]%
\hspace*{1em} \liuhao \linespread{0.8}\selectfont
\parbox{\textwidth}{
\heiti\makebox[\widthof{\makebox{*}实}][r]{实}验时间: #1\\
\heiti\makebox[\widthof{\makebox{*}实}][r]{报}告时间: #2\\
\heiti\makebox[\widthof{\makebox{*}实}][r]{$\dagger$指}导教师\\
\heiti\makebox[\widthof{\makebox{*}实}][r]{\makebox{*}学}号：#3 \\
\heiti\makebox[\widthof{\makebox{*}实}][r]{\makebox{*}E-}mail：#4 \\
}}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 首页页眉页脚定义
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\fancypagestyle{plain}{ \fancyhf{}\vspace{10pt}
\lhead{}
\chead{\vspace{0.6cm}\centering{{\heiti\xiaowuhao 清\ 华\ 大\ 学\ 近\ 代\ 物\ 理\ 实\ 验}\\[-0.5ex]
{{\xiaowuhao Experiments in Modern Physics, Tsinghua University}}}}
\rhead{}
\lfoot{}
\cfoot{}
\rfoot{}
\renewcommand{\headrule}{%
\hrule height0.4pt width \headwidth \vskip1.0pt%
\hrule height0.4pt width \headwidth \vskip-2pt}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\pagestyle{fancy} 
\fancyhf{}
\fancyhead[CE]{{\footnotesize 清华大学近代物理实验}}
\fancyhead[LE,RO]{~~~\footnotesize\thepage ~~~}
\lfoot{}
\cfoot{}
\rfoot{}


\DeclareMathSizes{10.5}{10}{6.8}{4.2}
\setlength{\abovedisplayskip}{2.5mm}
\setlength{\belowdisplayskip}{2.5mm}
