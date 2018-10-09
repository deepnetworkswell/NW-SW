clear all
% Input
%sequences to be aligned
Sequence1='QFTQEQQDIEDNCQC'             %Vertical
Sequence2='FFGQEQSYIEDNCND'         %Horizontal



%Scores
matchScore=1;            %match rewards
gapPenalty=0;           %gap penalty
mismatchPenalty=0;      %mismatch penalty


%determining the length of each sequence
lengthSequence1 = length(Sequence1);
lengthSequence2 = length(Sequence2);

%Initialize the first column
for i = 2:lengthSequence1+1
    nwScoringMatrix(i,1) = i*gapPenalty;
end
%Initialize the first row
for i = 2:lengthSequence2+1
    nwScoringMatrix(1,i) = i*gapPenalty;
end

for i = 2:lengthSequence1+1
    char1 = Sequence1(i-1);
    for j = 2:lengthSequence2+1
        char2 = Sequence2(j-1);
         if char1 == char2
            alignmentScore = matchScore;
         else 
            alignmentScore = mismatchPenalty;
         end
    
        score1 = nwScoringMatrix(i-1,j-1) + alignmentScore;       %diagonal
        score2 = nwScoringMatrix(i-1,j) + gapPenalty;             %up
        score3 = nwScoringMatrix(i,j-1) + gapPenalty;             %left
        scores=[score1 score2 score3];
        % two matrices to show the score and trDirection
        [maxScore,trDirection] = max(scores);
        nwScoringMatrix(i,j)=maxScore;
        switch trDirection
            case 1
                trDir='D';
            case 2
                trDir='U';
            case 3
                trDir='L';
        end
        nwTracebackMatrix(i-1,j-1)=trDir;
    end
end

% Creating sequences
alignedSequence1 = ''; alignedSequence2 = '';
i = lengthSequence1;
j = lengthSequence2;

% Tracebacking algorithm
while j > 0 && i > 0
       if nwTracebackMatrix(i,j)=='D'   
           alignedSequence1 =[Sequence1(i) alignedSequence1];
           alignedSequence2 =[Sequence2(j) alignedSequence2];
           i=i-1;
           j=j-1;
       elseif nwTracebackMatrix(i,j)=='U'   
           alignedSequence1 = [Sequence1(i) alignedSequence1];
           alignedSequence2 = ['-' alignedSequence2];
           i=i-1;
       elseif nwTracebackMatrix(i,j)=='L' 
           alignedSequence1 = ['-' alignedSequence1];
           alignedSequence2 = [Sequence2(j) alignedSequence2];
           j=j-1;
       end
end

% Output
score=nwScoringMatrix(lengthSequence1+1,lengthSequence2+1);
fprintf('Score = %i\n',score);
fprintf('Alignment = %s\n',alignedSequence2);
fprintf('Alignment = %s\n',alignedSequence1);