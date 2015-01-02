%nw aligment algorithm for HW1

%creates a matlab function
function [TotalScore, Alignment_A, Alignment_B] = nw_anchored(seq1, seq2, txt)
%function nw_anchored(seq1, seq2, txt)
%check number of parameters
%if (nargin == 2 || nargin == 3)
    %Take in inputs (seq1, seq2, [matches.txt])
    %[header, sequence]
    [Header1,A] = fastaread(seq1);
    %assign variables A & B to be the sequences
    [Header2,B] = fastaread(seq2);
    %check if there is a third parameter and load it
    if (nargin == 3)
        matches = load(txt);
    end    
%end
%Computing the matrix (F)
%gap penalty = -2
d = -2;
col = length(A);
row = length(B);
% rows|columns

%preallocating matrix F with the gap penalty as the top row and column
F = zeros(row+1,col+1);
F(2:end,1) = d * (1:row)';
F(1,2:end) = d * (1:col);

%match = 1, mismatch = -3
m = 1;
s = -3;

%variables for the wrapper function
if (nargin == 3)
    matchA = 1;
    matchB = 1;
    
    %sets the counters to be the first known match position
    counterA = matches(matchA);
    counterB = matches(matchB,3);
end

TotalScore = 0;
%scores = zeros(row, col);
%Filling in the matrix
for i=2:row+1
    for j=2:col+1        
        %wrapper function for input matches file
        if (nargin == 3)
           %j is columnA for human(needs to be input first)
           if (j == counterA && i == counterB)
               %assigns the position as a match
               F(i,j) = F(i-1,j-1) + scores(A(j-1),B(i-1));
               stop = matches(matchA,2);
               final = size(matches,1);
               %checks if we've reached the end of the known match
               if (counterA == stop)
                   %checks if we've reached the end of the text file
                   if (matchA == final)
                       continue
                   else
                       matchA = matchA + 1;
                   end
               else
                   %counts up through the known matched sequence
                   counterA = counterA + 1;
                   counterB = counterB + 1;
               end
               %skips to the next iteration of the for loop
               continue
            end
        end
        %j is going through the columns which is A
        %if the two positions match, the index in scores is listed as such
        if (A(j-1) == B(i-1))
            scores(A(j-1),B(i-1)) = m;
        else
            scores(A(j-1),B(i-1)) = s;
        end
        
        %Filling-in partial alignments here  
        Match     = F(i-1,j-1) + scores(A(j-1),B(i-1));
        MismatchA = F(i, j-1) + d;
        MismatchB = F(i-1, j) + d;
        %computing the final score of the alignment and assigning it to F
        Temp = [Match MismatchA MismatchB];
        F(i,j) = max(Temp);
          % end
    end
end
%{
%matched scores
if (nargin == 3)
    %counts through number of rows in the text file
    for q=1:size(matches,1)
        %sets the start and end points of the human sequence
        start_1 = matches(q);
        end_1 = matches(q,2);
        %counts the total number of matches and adds it to total score
        total = (end_1 - start_1) + 1;
        TotalScore = TotalScore + total;
    end
end
%}
%}
%traceback part
Alignment_A = '';
Alignment_B = '';
i = length(B)+1; %row
j = length(A)+1; %col

while (i>1 && j>1)
   Score = F(i,j);
   DIAG = F(i-1,j-1);
   %LEFT = F(i-1,j);
   UP   = F(i,j-1);
   
   %if scores are equal to the diagonal, there is no gap.
   if (Score == DIAG + scores(A(j-1),B(i-1)))
       Alignment_A = strcat(Alignment_A, A(j-1));
       Alignment_B = strcat(Alignment_B, B(i-1));
       %computes score, checks to see if the alignment are the same
       %characters
       if (A(j-1) == B(i-1))
           TotalScore = TotalScore + 1;
       else
           %mismatch
           TotalScore = TotalScore + s;
       end
       i = i-1;
       j = j-1;
   
   %gap in sequence B
   elseif (Score == UP + d)
       Alignment_A = strcat(Alignment_A, A(j-1));
       Alignment_B = strcat(Alignment_B, '-');
       j = j-1;
       %computes score by adding the gap penalty
       TotalScore = TotalScore + d;
   %gap in sequence A
   else
       Alignment_A = strcat(Alignment_A, '-');
       Alignment_B = strcat(Alignment_B, B(i-1));
       i = i-1;
       TotalScore = TotalScore + d;
   end
end
%If at the end of one sequence fills in the rest with gaps
while(j>1)
   Alignment_A = strcat(Alignment_A, A(j-1));
   Alignment_B = strcat(Alignment_B, '-');
   j = j-1;
   %-2 for gaps
   TotalScore = TotalScore - 2;
end
while(i>1)
   Alignment_A = strcat(Alignment_A, '-');
   Alignment_B = strcat(Alignment_B, B(i-1));
   i = i-1;
   TotalScore = TotalScore - 2;
end

%displays the alignment score and alignments
disp('Alignment Score =');
disp(TotalScore);
disp(Alignment_A);
disp(Alignment_B);
