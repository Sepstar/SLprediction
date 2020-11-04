[header,sequence] = fastaread('protein.fasta');
i = 1;
score = zeros(1250,1250);
for i = 1:1250
    for j = (i+1):1250
        score(i,j)=swalign(sequence{i},sequence{j})/(sqrt(swalign(sequence{i},sequence{i}))*sqrt(swalign(sequence{j},sequence{j})));
    end
    i
end
%xlswrite('SLscore.xlsx',score);
scoreall = score + score';
for i=1:length(scoreall)
    for j=1:length(scoreall)
        if i==j
            scoreall(i,j)=1;
        end
    end
end
xlswrite('SLscoreall.xlsx',scoreall);
xlswrite('SLheader.xlsx',header);
