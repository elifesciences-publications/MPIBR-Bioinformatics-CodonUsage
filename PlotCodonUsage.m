%% PlotCodonUsage
clc
clear variables
close all

%% read codon table
data = readCodonTable('codonTable_mouse_05Aug2017.txt');

%% read expression weights
wout = readExpressionWeights('mouseHippocampus_expressionWeights.txt');

%% read query list
list = readQueryList('list_kinase.txt');

%% intersect weights with transcripts
[symunq, ~, idx_unq] = unique(data.symbols);
[~, idx_fwd, idx_rev] = intersect(symunq, wout.symbols);
ww = zeros(max(idx_unq), 1);
ww(idx_fwd) = wout.weights(idx_rev);
ww = ww(idx_unq);

%% intersect kinases with transcripts
[symunq, ~, idx_unq] = unique(data.symbols);
[~, idx_fwd] = intersect(symunq, list.symbols);
idx_kinases = false(length(symunq), 1);
idx_kinases(idx_fwd) = true;
idx_kinases = idx_kinases(idx_unq);

%% calculate count table
counts = bsxfun(@times, data.counts, ww);

counts_total = sum(counts, 1);
w = zeros(length(counts_total), 1);
[~, ~, idx_group] = unique(data.lists.AA);
for g = 1 : max(idx_group)
    
    %% calculate for the total
    w_aa = counts_total(idx_group == g) ./ sum(counts_total(idx_group == g));
    w(idx_group == g) = w_aa ./ max(w_aa);
    
end

cai = exp(sum(bsxfun(@times, counts, log(w')), 2) ./ (sum(counts, 2) - 1));


%{
idx_base = ww > 0;
idx_query = strcmp(data.symbols, 'Eif2ak1');

[f_all, x_all] = ecdf(cai(idx_base));
[f_kin, x_kin] = ecdf(cai(idx_kinases & idx_base));

figure('color','w');
hold on;
plot(x_all, f_all, 'k', 'linewidth', 1.2);
plot(x_kin, f_kin, 'color', [30,144,255]./255, 'linewidth', 1.2);
hold off;

%
[N_all, E_all] = histcounts(cai(idx_base), 100);
[N_kin, E_kin] = histcounts(cai(idx_base & idx_kinases));


cai_query = mean(cai(idx_query));
cai_median = median(cai(idx_base));
cai_kinases = median(cai(idx_base & idx_kinases));

% inverse percentile
X = cai(idx_base);
Q = cai_query;
nless = sum(X < Q);
nequal = sum(X == Q);
centile = 100 * (nless + 0.5*nequal) / length(X);
%}

%{
figure('color','w');
hold on;
h(1) = plot(E_all(1:end-1), N_all./max(N_all),'k','linewidth',1.2);
%plot([cai_median, cai_median], [0, 1], 'k');
h(2) = plot(E_kin(1:end-1), N_kin./max(N_kin),'color',[30,144,255]./255,'linewidth',1.2);
%plot([cai_kinases, cai_kinases], [0, 1], 'color',[30,144,255]./255);
h(3) = plot([cai_query, cai_query],[0,1],'r','linewidth',1.2);
hold off;
hl = legend(h, sprintf('genes %d', length(unique(data.symbols(idx_base)))),...
               sprintf('kinases %d',length(unique(data.symbols(idx_base & idx_kinases)))),...
               'Eif2ak1');
set(hl,'edgecolor','w','location','northwest');
title('Genes expressed in hippocampal slices','fontweight','normal','fontsize',14);
xlabel('codon adaptation index','fontsize',12);
ylabel('relative count', 'fontsize',12);
print(gcf,'-dpng','-r300','figureX_CAI_Eif2ak1.png');

%}








%% --- FUNCTIONS --- %%%
function out = readQueryList(file_name)
    fh = fopen(file_name, 'r');
    txt = textscan(fh, '%s %s', 'delimiter', '\t');
    fclose(fh);
    out.symbols = txt{1};
    out.name = txt{2};
end


function out = readExpressionWeights(file_name)
    fh = fopen(file_name, 'r');
    txt = textscan(fh, '%s %n', 'delimiter', '\t');
    fclose(fh);
    out.symbols = txt{1};
    out.weights = (txt{2}.*1e6)./sum(txt{2});
    idx_filter = out.weights < 2;
    out.symbols(idx_filter) = [];
    out.weights(idx_filter) = [];
end


function data = readCodonTable(file_name)

    fmt_header = repmat({'%s'}, 66, 1);
    fmt_header(1:2) = {'%*s'};
    fmt_header = sprintf('%s ', fmt_header{:});
    fmt_header(end) = [];

    fmt_counts = repmat({'%n'}, 66, 1);
    fmt_counts(1:2) = {'%s'};
    fmt_counts = sprintf('%s ', fmt_counts{:});
    fmt_counts(end) = [];

    fh = fopen(file_name, 'r');
    list_codons = textscan(fh, fmt_header, 1, 'delimiter', '\t');
    list_aa = textscan(fh, fmt_header, 1, 'delimiter', '\t');
    list_abr = textscan(fh, fmt_header, 1, 'delimiter', '\t');
    txt = textscan(fh, fmt_counts, 'delimiter', '\t');
    fclose(fh);
    data.lists.codons = [list_codons{:}]';
    data.lists.AA = [list_aa{:}]';
    data.lists.name = [list_abr{:}]';
    data.symbols = txt{1};
    data.transcripts = txt{2};
    data.counts = [txt{3:end}];

end





