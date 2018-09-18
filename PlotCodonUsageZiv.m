%% PlotCodonUsageZiv
clc
clear variables
close all


%% read codon table
data = readCodonTable('codonTable_mouse_05Aug2017.txt');

%% read expression weights
wout = readExpressionWeights('mouseHippocampus_expressionWeights.txt');
%% intersect weights with transcripts
[symunq, ~, idx_unq] = unique(data.symbols);
[~, idx_fwd, idx_rev] = intersect(symunq, wout.symbols);
ww = zeros(max(idx_unq), 1);
ww(idx_fwd) = wout.weights(idx_rev);
ww = ww(idx_unq);


%% read query list
list = readQueryList('EV1_ZivPaper2013/list_EV1.txt');


%% intersect list with transcripts
[symunq, ~, idx_unq] = unique(data.symbols);
[~, idx_fwd, idx_rev] = intersect(symunq, list.symbols);
idxZiv = false(length(symunq), 1);
idxZivQuery = false(length(symunq), 1);
idxZivSynaptic = false(length(symunq), 1);
idxZiv(idx_fwd) = true;
idxZivQuery(idx_fwd) = list.query(idx_rev);
idxZivSynaptic(idx_fwd) = list.synaptic(idx_rev);
idxZiv = idxZiv(idx_unq);
idxZivQuery = idxZivQuery(idx_unq);
idxZivSynaptic = idxZivSynaptic(idx_unq);

%% calculate CAI
ww = ones(length(data.counts),1);
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



%% test
%{
V = [cai;...
     cai(idxZiv);...
     cai(idxZivSynaptic);...
     cai(idxZivQuery)];
G = [ones(length(cai),1);...
     2*ones(sum(idxZiv),1);...
     3*ones(sum(idxZivSynaptic),1);...
     4*ones(sum(idxZivQuery),1)];
 
[h,p,s] = kruskalwallis(V, G);
multcompare(s)
%}

%% plot results
[fAll,xAll] = hist(cai, 100);
[fZiv,xZiv] = hist(cai(idxZiv), 100);
[fSyn,xSyn] = hist(cai(idxZivSynaptic), 100);
[fQry,xQry] = hist(cai(idxZivQuery), 100);

figure('color','w');

h(1) = subplot(2, 2, 1);
hb = bar(xAll, fAll, 'facecolor', [.45, .45, .45],...
                'edgecolor', [.45, .45, .45],...
                'facealpha', 0.3);
hl = legend(hb,'NCBI annotated proteome'); 
set(hl,'edgecolor','w','location','northoutside');
xlabel('codon adaptation index');
ylabel('number of proteins');
            
h(2) = subplot(2, 2, 2);
hb = bar(xZiv, fZiv, 'facecolor', [30, 144, 255] ./ 255,...
                'edgecolor', [30, 144, 255] ./ 255,...
                'facealpha', 0.3);
hl = legend(hb,'Hakim et.al. 2013 : proteome'); 
set(hl,'edgecolor','w','location','northoutside');      
xlabel('codon adaptation index');
ylabel('number of proteins');            
            
h(3) = subplot(2, 2, 3);
hb = bar(xSyn, fSyn, 'facecolor', [50, 205, 50] ./ 255,...
                'edgecolor', [50, 205, 50] ./ 255,...
                'facealpha', 0.3);
hl = legend(hb,'Hakim et.al. 2013 : synaptic'); 
set(hl,'edgecolor','w','location','northoutside');            
xlabel('codon adaptation index');
ylabel('number of proteins');            
            
h(4) = subplot(2, 2, 4);
hb = bar(xQry, fQry, 'facecolor', [148,0,211] ./ 255,...
                'edgecolor', [148,0,211] ./ 255,...
                'facealpha', 0.3);
hl = legend(hb,'Hakim et.al. 2013 : paradoxical synthesis'); 
set(hl,'edgecolor','w','location','northoutside');
xlabel('codon adaptation index');
ylabel('number of proteins');            
            
set(h,'box', 'off',...
      'xlim', [0.6, 1]);

print(gcf, '-dpng', '-r300', 'figure_Hakim2013_CAI.png');




%% functions
function out = readQueryList(file_name)
    fh = fopen(file_name, 'r');
    txt = textscan(fh, '%s %n %n', 'delimiter', '\t');
    fclose(fh);
    out.symbols = txt{1};
    out.query = logical(txt{2});
    out.synaptic = logical(txt{3});
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

