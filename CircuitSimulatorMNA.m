%Here we prompt user for the file name.
%nodevoltages function takes a file name in txt form as an input,extension
%of the file must be indicated. Returns an array of node voltages and
%prints out the voltage values of each node.

file = input('Enter the name of the file: ','s');
nodes = nodevoltages(file);

function nodes = nodevoltages(file)
%%nodevoltages function takes a file name in txt form as an input,extension
%%of the file must be indicated. Returns an array of node voltages and
%%prints out the voltage values of each node.
%%see also modifiednodevoltages

%Reading the file
fid= fopen(file,'r');
FormatSpec = '%s %d %d %d';
%We will put the read cell into CircuitMat cell array
CircuitMat = textscan(fid,FormatSpec);
fclose(fid);
%We have put our text file columns into a cell.
%% 
%We must group the elements in order to find the number of elements of
%each type
CircuitMat{1,4} = double(CircuitMat{1,4}); %We change the integer values to double precision values.
[voltagesources,resistances,currentsources] = group(CircuitMat);%group function groups the elements and returns 
%the arrays of the elements

%we find unique values of nodes, by adding the column vectors into a
%greater column so that we allow unique function to return the unique values of
%both columns designated for nodes in CircuitMat, and with size function we
%get its size we subtract one because 0 node is not included in node voltages 
%hence number of nodes.
%n is the number of nodes
n = size(unique([CircuitMat{2};CircuitMat{3}]),1)-1;
%m is the number of voltage sources
m = size(voltagesources,1);
%r is the number of resistances
r = size(resistances,1);
%inum is the number of current sources
inum= size(currentsources,1);
%%
%We will create 4 sub matrices G B C D to create the matrix A (m+n)x(m+n)
%n is the number of nodes and m is the number of independent voltage
%sources

%%
%G matrix is created like whenever there is a resistor connected to node i
%the conductance of resistor (1/R) will appear ith row ith column to be
%summed and whenever there are resistors between nodes, lets say node i and
%j, their counductance will appear in Gij and Gji as negative.

G = zeros(n,n);
%with this for loop we fill all entries, including diagonal ones but
%diagonal ones will be replaced by the for loop below.
for i = 1:r
    if resistances{i,2} ~= 0 && resistances{i,3} ~= 0 %resistances which are conencted to 0 node wont
        %appear in the non-diagonal entries 
        G(resistances{i,2},resistances{i,3}) = G(resistances{i,2},resistances{i,3})-1/resistances{i,4};
        G(resistances{i,3},resistances{i,2}) = G(resistances{i,3},resistances{i,2})-1/resistances{i,4};
    end
end
%with this for loop we managed to fill the diagonals of G matrix with
%correct values.
for i = 1:n 
    summ = 0;
    %for loop below checks if a resistance is conencted to a node other
    %than 0 it will appear in the diagonal, i is the node number that is
    %being checked.
    for j = 1:r % r is the number of resistances
        if resistances{j,2} == i 
            summ = summ + 1.0/resistances{j,4};
        elseif resistances{j,3} == i
            summ = summ + 1.0/resistances{j,4};
        end
    end
    G(i,i) = summ;
end


%%
%Now we will create B matrix, it is an mxn matrix, columns correspond to
%voltage sources and rows correspond to nodes, if the positive terminal
%of the voltage source i is connected to node k then the element B(k,i) =
%1, if negative terminal is connected to k B(k,i) = -1 otherwise it is
%zero.
B = zeros(n,m);
%positive terminals of current and voltage sources are on the third column
%of source cell arrays.
for i = 1:m
    if voltagesources{i,2}~= 0
        B(voltagesources{i,2},i) = -1;
    end
    if voltagesources{i,3}~= 0
        B(voltagesources{i,3},i) =1;
    end
end
%%
%C matrix is just the transpose of the matris B since there is no dependent
%sources
C = transpose(B);
%%
%D matrix is a mxm matrix consists of zeros entirely
D = zeros(m,m);
%%
A = [ G B ; C D];
% we have built our matrix A
%%
% Now we will build our matrix x, x hold the unknown quantities, and is
% built by 2 smaller matrice v and j. x is a (m+n)x1 matrix
%%
%v matrix is nx1 matrix holding unknown node voltages
v = sym('v_',[n 1]);
%%
%j matrix is mx1 matrix holding unknown currents passing through voltage
%sources.
J = sym('I_V',[m 1]);
%%
%now we constitute the matricex v and J into x matrix
x = [ v;J];
%%
%we will now create the z matrix,
%z matrix (m+n)x1 and it will be constituted by 2 smaller matrices I and e 
%%
%I matrix will be nx1 matrix of entering currents of
%independent current sources to a certain node.
I = zeros(n,1);
for i= 1:inum
    if currentsources{i,3} ~= 0
        I(currentsources{i,3},1) = currentsources{i,4};
    end
end
%%
% e matrix will be mx1 matrix holding the voltage values of voltage sources
e = zeros(m,1);
for i = 1:m
    e(i,1) = voltagesources{i,4};
end
%%
% Now we will constitute I and e matrices into z matrix
z = [I;e];
%%
%Now we have all the matrices needed we can simply solve these linear
%equations system.
 x = linsolve(A,z);
%%
%Now we will create the output of our code.
nodes = [];
for i = 1:n
    nodes(end+1)=x(i);
end
for i = 1:n
    fprintf('The voltage of node %d is %.4fV \n',i,x(i));
end
end
function [voltagesources,resistances,currentsources] = group(cell)
%In this function we use strncmp to compare the first letter of the
%elements to identify which group they belong.
%the outputs will be cellarrays voltagesources,resistances,currentsources and input will be a cell

%this for loop below will create a cell only involving resistances
j = 1;
for i = 1:size(cell{1})
    if strncmp('R',cell{1}(i),1)
        resistances{j,1} = cell{1}(i);
        resistances{j,2} = cell{2}(i);
        resistances{j,3} = cell{3}(i);
        resistances{j,4} = cell{4}(i);
        j = j+1;
    end
end
if j ==1
    resistances = {};
end
%this loop below will group Voltage Sources
j =1;
for i = 1:size(cell{1})
    if strncmp('V',cell{1}{i},1)
        voltagesources{j,1} = cell{1}(i);
        voltagesources{j,2} = cell{2}(i);
        voltagesources{j,3} = cell{3}(i);
        voltagesources{j,4} = cell{4}(i);
        j = j+1;
    end
end
if j ==1
    voltagesources = {};
end
%this loop will group current sources
j = 1;
for i = 1:size(cell{1})
    if strncmp('I',cell{1}{i},1)
        currentsources{j,1} = cell{1}(i);
        currentsources{j,2} = cell{2}(i);
        currentsources{j,3} = cell{3}(i);
        currentsources{j,4} = cell{4}(i);
        j = j+1;
    end
end
if j==1
    currentsources = {};
end
end