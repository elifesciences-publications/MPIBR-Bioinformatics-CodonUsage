%% PlotCodonUsage
clc
clear variables
close all

%% read codon table
fmt_header = repmat({'%s'}, 66, 1);
fmt_header(1:2) = {'%*s'};
fmt_header = sprintf('%s ', fmt_header{:});
fmt_header(end) = [];

fmt_counts = repmat({'%n'}, 66, 1);
fmt_counts(1:2) = {'%s'};
fmt_counts = sprintf('%s ', fmt_counts{:});
fmt_counts(end) = [];

fh = fopen('codonTable_mouse_05Aug2017.txt', 'r');
list_codons = textscan(fh, fmt_header, 1, 'delimiter', '\t');
list_codons = [list_codons{:}]';
list_aa = textscan(fh, fmt_header, 1, 'delimiter', '\t');
list_aa = [list_aa{:}]';
list_abr = textscan(fh, fmt_header, 1, 'delimiter', '\t');
list_abr = [list_abr{:}]';
txt = textscan(fh, fmt_counts, 'delimiter', '\t');
fclose(fh);
symbols = txt{1};
transcripts = txt{2};
counts = [txt{3:end}];

