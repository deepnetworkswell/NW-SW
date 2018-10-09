clear all
% Input
%sequences to be aligned
Sequence1= 'MTTVVTPNPGLGKAS'    %Vertical
Sequence2= 'VSTVVLENPGLGRAL'      %Horizontal
%Scores
matchScore=1;           %match rewards
gapPenalty=1/3;           %gap penalty
mismatchPenalty=1/3;      %mismatch penalty
    
lengthSequence1 = length(Sequence1);
lengthSequence2 = length(Sequence2);


%Initializing first column and first row
for n = 2:lengthSequence1+1
    swScoringMatrix(n,1) = 0;
end

for n = 2:lengthSequence2+1
    swScoringMatrix(1,n) = 0;
end

%Diagonally
for i = 2:lengthSequence1+1
    char1 = Sequence1(i-1);
    for j = 2:lengthSequence2+1
        char2 = Sequence2(j-1);
         if char1 == char2
            alignmentScore = matchScore;
        else 
            alignmentScore = -mismatchPenalty;
        end
        score1=0;
        score2 = swScoringMatrix(i-1,j-1) + alignmentScore;       %diagonal
        scores=[score1 score2];
        % Matrix to record the score and direction
        [maxScore,trDirection] = max(scores);
        swScoringMatrix(i,j)=maxScore;
        swTracebackMatrix(i-1,j-1)=trDirection;
    end
end


%Horizontally
for i = 2:lengthSequence1+1
    for j = 2:lengthSequence2+1
        oldValueij=swScoringMatrix(i,j);
        hscore=oldValueij-1;
        k=j+1;
        while (hscore>0 && k<=lengthSequence2+1)
            hscore=hscore-0.333333;
            if(swScoringMatrix(i,k)<hscore)
                swScoringMatrix(i,k)=hscore;
                swTracebackMatrix(i-1,k-1)=3;
            end
            k=k+1;    
        end
    end
end

    
 %Vertically
for j = 2:lengthSequence2+1
    for i = 2:lengthSequence1+1

        oldValueij=swScoringMatrix(i,j);
        hscore=oldValueij-1;
        k=i+1;
        while (hscore>0 && k<=lengthSequence1+1)
            hscore=hscore-0.333333;
            if(swScoringMatrix(k,j)<hscore)
                swScoringMatrix(k,j)=hscore;
                swTracebackMatrix(k-1,j-1)=4;
            end
            k=k+1;    
        end
    end
end
    

% Traceback
[M,I] = max(swScoringMatrix(:));
[I_row, I_col] = ind2sub(size(swScoringMatrix),I)

% Sequences are built
alignedSequence1 = '';  alignedSequence2 = '';
i = I_row-1;
j = I_col-1;


while j > 0 && i > 0

       % Diagonal
       if swTracebackMatrix(i,j)==2   
           alignedSequence1 =[Sequence1(i) alignedSequence1];
           alignedSequence2 =[Sequence2(j) alignedSequence2];
           i=i-1;
           j=j-1;
       % Up
       elseif swTracebackMatrix(i,j)==4  
           alignedSequence1 = [Sequence1(i) alignedSequence1];
           alignedSequence2 = ['-' alignedSequence2];
           i=i-1;
       % left
       elseif swTracebackMatrix(i,j)==3   
           alignedSequence1 = ['-' alignedSequence1];
           alignedSequence2 = [Sequence2(j) alignedSequence2];
           j=j-1;
       % Zero
       elseif swTracebackMatrix(i,j)==1  
           break;
       end    
end
% Output
[score,I] = max(swScoringMatrix(:));
fprintf('Score = %i\n',score);
fprintf('Alignment = %s\n',alignedSequence2);
fprintf('Alignment = %s\n',alignedSequence1);



