function y = fano_algorithm()
%g1 is the first generator function which can be hardcoded according to
%preference.
%g2 is the second generator function which can be hardcoded according to
%preference of generator function.
g1 = [1 1 1 0 1 0];
g2 = [1 1 1 1 0 1];
%A array stores the matrix contaning bit shifted values for the given
%encoders for multiplication with the dataword to get the codeword
A = zeros(8,26); 
%x is the transition table to store all the values of possible transition
%for the given encoder for input 0 and 1.The first 5 bits of the transtion
%table contains all the possible commbination of 5 bit current states for
%the code.The 6th bit is the current input which can be either 0 or 1.The
%bits from 7 to 11 are the bits corresponding to the next state transition
%from the current input state of bits. The last two bits that is 12 and 13
%stands for the output given by using the two generator functions.
x = zeros(64,13);

%this part of the code is responsible for encoding the provided data word
%-------------Encoding----------------------------------------------------
k=1;
for i=1:6
  p(k) = g1(i);
  p(k+1) = g2(i);
  k=k+2;
end

k=1;
for i=1:8,
   A(i,[k:k+11]) = p;
   k=k+2;
end;
fprintf('\nMENU\n');
inputdata = 'Enter the data word (8 bits): ';
data = input(inputdata,'s');
data = data - '0';
encoded = zeros(1,26);
encoded = mod(data*A,2);
%-------------Encoding part ends here-------------------------------------
%now we move on to the generation of the transition table for enoding
%--------------------Generation of the transition table------------------
k=1;
for i=0:31
   int_array(k) = i;
   int_array(k+1) = i;
   k=k+2;
end
%converting that to a binary array containing only binary numbers 
bin_array = dec2bin(int_array,5) - '0';
%creating all possible combinations of the binary number in the transition table
x(:,[1:5]) = bin_array;
%this is for filling the input column of the matrix
alt = 0;
%this for loop fills all the bits of the 1x64 transition table 
for i=1:64
   for j=8:11
      x(i,j) = x(i,j-7);
   end
   
   x(i,6) = alt;
   if alt==1
      alt = 0;
   else
      alt = 1;
   end
   x(i,7) = x(i,6);
   temp1 = x(i,[1:5]);
   temp2 = [x(i,6),temp1];
   x(i,12) = mod(sum(g1.*temp2),2);
   x(i,13) = mod(sum(g2.*temp2),2);
end
fprintf('\n1. Enter received codeword and then the number of error bits between corrected codeword and actual sent codeword is displayed\n2. Display a bar graph that displays the error correction and detection percentages\n');
choice = input('\nEnter the choice: '); 
switch choice
    case 2    
%-------------end of creation of the transition table----------------------
% This part deals with the decoding using fano algorithm for all thresholds
% from 2 to specified global range.
disp('Processing.........');
%inarr is the array to which we introduce errors to the encoded word before
%passing it to the fano algorithm
inarr = encoded;
%max_errors is the global variable which contains the limit to the number
%of errors introduced in the code word by the system.
max_errors = 4;
%threshold_limit contains the limit upto which the subplots will be
%plotted.
threshold_limit = 5;
%index variable is used to traverse the input array or inarr 
index = 1;
%outarr contains the sequence of output bits found out by the fano
%algorithm.
outarr = zeros(1,26);
m = 1;
%error_detection is an array that stores the no of errors detected where
%the no of errors detected for i bit errors introduced in the codeword is
%stored in location indexed by i in the array.
error_detection = zeros(1,26);
%error_correction is an array that stores the no of errors corrected where
%the no of errors corrected for i bit errors introduced in the codeword is
%stored in location indexed by i in the array.
error_correction = zeros(1,26);
%error_detection_avg is responsible for calculating the average number of
%errors detected
error_detection_avg = zeros(1,26);
%error_correction_avg finds the average number of errors corrected among
%the combinations generated
error_correction_avg = zeros(1,26);
%comb array contains all the indexes that range from 1 to 26 of indices so
%that nchoosek can chose combination among these
comb = [1:26];
for th = 2:threshold_limit
%note: the i varible here stands for the number of errors introduced into
%the codeword by the system/program,here we restrict to maximum of 4 errors
%introduced due to timeconstraint and this value can be changed globally.
for i = 1:max_errors
    temp = nchoosek(comb,i);
    [r,c] = size(temp);
    for j = 1:r
        for k = 1:c
            inarr(1,temp(j,k)) = ~inarr(1,temp(j,k));
        end
        [check,outarr] = trav(m,x,index,inarr,outarr,encoded,th);
        if check == 0
            error_detection(1,i) = error_detection(1,i) + 1;
        elseif sum(xor(outarr,encoded)) == 0
            error_detection(1,i) = error_detection(1,i) + 1;
            error_correction(1,i) = error_correction(1,i) + 1;
        elseif sum(xor(outarr,encoded)) ~= 0
            if sum(xor(inarr,outarr)) ~= 0  %False detection
               error_detection(1,i) = error_detection(1,i) + 1;
            end
        end
        inarr = encoded;
        index = 1;
        m = 1;
       
    end 
    error_detection_avg(1,i) = (error_detection(1,i)/r)*100;
    error_correction_avg(1,i) = (error_correction(1,i)/r)*100;
end
%error_detection
%error_correction
%error_detection_avg
%error_correction_avg

%the graphin array is responsible for creating a 4row 2 column array
%containing the four rows as pairs of 4 error correction and detection for
%errors introduced upto 4 bits due to time constraint
graphin = zeros(4,2);
for i = 1:4
    graphin(i,1) = error_detection_avg(1,i);
    graphin(i,2) = error_correction_avg(1,i);
end
%note: figure which is commented below can be used without the sub-plot
%statement to plot seperate graphs for all the threshold values specified
%upto the limit of theshold values.
%figure;

%the subplot statement here plots threshold_limit -1 number of graphs as it
%starts from 2
subplot(1,threshold_limit-1,th-1);
%this statement plots the bargraph
bar(graphin)
ylim([0 130])
title(['Convolution Coding(Threshold=' num2str(th) ')'])
legend('Error detection (%)','Error correction (%)')
xlabel('Number of errors')
ylabel('Percentage (%)')
y = 0;
error_detection = zeros(1,26);
error_correction = zeros(1,26);

end
    case 1 
        encoded
%encoded array holds encoded codeword.          
        threshold = input('Enter the threshold: ');
%User enters the threshold value        
        outarr = zeros(1,26);
        inarr = encoded;
%inarr is the array holding the received codeword
        posno = input('Enter the number of positions where bits have to be complemented: ');
        for i = 1:posno
            pos = input(['Enter position ' num2str(i) ': ']);
            inarr(1,pos) = ~inarr(1,pos);
        end
%Errors are introduced at user entered positions in inarr array
        [check,outarr] = trav(1,x,1,inarr,outarr,encoded,threshold);
        if check == 0
            fprintf('Error is detected\nError is not corrected\nCorrected codeword is not generated by this algorithm\n');
        elseif sum(xor(outarr,encoded)) == 0
            fprintf('\nError is detected and corrected\n');
            outarr
            error = floor(check/10)
        elseif sum(xor(outarr,encoded)) ~= 0
            if sum(xor(inarr,outarr)) ~= 0
                fprintf('\nError is detected but not corrected\n');
                outarr
                error = floor(check/10)
            else
                fprintf('The algorithm has falsely declared the input to be correct because it corresponds\nto another valid codeword\n');
                fprintf('Error is neither detected nor corrected\n');
                error = floor(check/10)
            end
        end
end          
end
%--------------------The error correction and detection analysis part ends
%here------------------------------------------------------------------


function [y,outarr] = trav(j,x,index,inarr,outarr,encoded,threshold)

%trav is a recursive function that implements sequential decoding algorithm
%to correct the received codeword. 
%Arguments of trav function : 
%1. x = State transition table
%2. index = Used to mark current state in transition table (determined by
%row number)
%3. inarr = Erroneous received codeword
%4. outarr = Corrected codeword by sequential decoding algorithm
%5. encoded = The actual encoded codeword of the dataword
%6. threshold = When the trav function is invoked by the calling program
%this variable is set to the threshold value. Every time a wrong path is 
%taken, the value of this variable decrements by 1. 

if threshold < 0
y = 0;
return
end

%If threshold value goes below 0, this means that the the function has 
%exceeded the threshold limit. Thus, we backtrack one step and return 0 to
%the calling function, which means that no solution has been reached from
%this path.

if threshold >= 0 && j < 26

%If the threshold limit has not yet been reached and the received codeword
%is not fully read, the statements following this if condition will be
%executed.
first =or(xor(x(index,12),inarr(1,j)),xor(x(index,13),inarr(1,j+1)));
second = or(xor(x(index+1,12),inarr(1,j)),xor(x(index+1,13),inarr(1,j+1)));

%first = 0 , if the present 2 bits from the received codeword under 
%consideration are same as the output bits from the present state of the
%state transition table with input bit 0.
%Otherwise, first = 1.
%first is the first outgoing branch.
%second = 0 , if the present 2 bits from the received codeword under 
%consideration are same as the output bits from the present state of the
%state transition table with input bit 1.
%Otherwise, second = 1.
%second is the second outgoing branch. 
if first == 0
%The next state from present state when input bit is 0 is taken.
   outarr(1,j:j+1) = x(index,12:13);
   newindex = power(2,4)*x(index,7) + power(2,3)*x(index,8) + power(2,2)*x(index,9) + power(2,1)*x(index,10) + power(2,0)*x(index,11);
%2*newindex + 1 = Row index of the next state in state transition table.   
   [y,outarr] = trav(j+2,x,2*newindex+1,inarr,outarr,encoded,threshold);
   if(y >= 1)
%If solution has been found, return to the previous calling function
     return
   else
%If solution is not found, go to the second branch.
     outarr(1,j:j+1) = x(index+1,12:13);
     newindex = power(2,4)*x(index+1,7) + power(2,3)*x(index+1,8) + power(2,2)*x(index+1,9) + power(2,1)*x(index+1,10) + power(2,0)*x(index+1,11);
%Since a wrong path is taken, the threshold variable is decremented by 1.     
     [y,outarr] = trav(j+2,x,2*newindex+1,inarr,outarr,encoded,threshold-1);
     return
   end
end

if second == 0   
%The next state from present state when input bit is 1 is taken.
   outarr(1,j:j+1) = x(index+1,12:13);
   newindex = power(2,4)*x(index+1,7) + power(2,3)*x(index+1,8) + power(2,2)*x(index+1,9) + power(2,1)*x(index+1,10) + power(2,0)*x(index+1,11);
%2*newindex + 1 = Row index of the next state in state transition table.
   [y,outarr] = trav(j+2,x,2*newindex+1,inarr,outarr,encoded,threshold);
   if(y >= 1)
%If solution has been found, return to the previous calling function
     return
   else
%If solution is not found, go to the first branch.
     outarr(1,j:j+1) = x(index,12:13);
     newindex = power(2,4)*x(index,7) + power(2,3)*x(index,8) + power(2,2)*x(index,9) + power(2,1)*x(index,10) + power(2,0)*x(index,11);
%Since a wrong path is taken, the threshold variable is decremented by 1.
     [y,outarr] = trav(j+2,x,2*newindex+1,inarr,outarr,encoded,threshold-1);
     return
   end
end

if first * second == 1   
%Both branches do not match the present two bits under consideration from
%the received codeword.
   outarr(1,j:j+1) = x(index,12:13);
   newindex = power(2,4)*x(index,7) + power(2,3)*x(index,8) + power(2,2)*x(index,9) + power(2,1)*x(index,10) + power(2,0)*x(index,11);
%First branch is taken first.
%2*newindex + 1 = Row index of the next state in state transition table.
%Since a wrong path is taken, the threshold variable is decremented by 1.
   [y,outarr] = trav(j+2,x,2*newindex+1,inarr,outarr,encoded,threshold-1);
   if(y >= 1)
%If solution has been found, return to the previous calling function.
     return
   else
%If solution is not found, go to the second branch.
     outarr(1,j:j+1) = x(index+1,12:13);
     newindex = power(2,4)*x(index+1,7) + power(2,3)*x(index+1,8) + power(2,2)*x(index+1,9) + power(2,1)*x(index+1,10) + power(2,0)*x(index+1,11);
%Since a wrong path is taken, the threshold variable is decremented by 1.
     [y,outarr] = trav(j+2,x,2*newindex+1,inarr,outarr,encoded,threshold-1);
     return
   end
end
    
end
%The following code will be executed only when the previous two if
%conditions evaluate to false. This means that threshold limit has not been
%reached and the entire received codeword is read. 
y = 1;
error = sum(xor(outarr,encoded));
%error = Hamming distance between actual codeword and the corrected
%codeword.
y = y + 10 * error;
%The ones place of y is 1 which signifies success. The error is obtained by
%truncating the last digit of y.
end