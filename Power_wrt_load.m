file = input('Enter the name of the file: ','s');
fid= fopen(file,'r');
FormatSpec = '%s %d %d %d';
CircuitMat = textscan(fid,FormatSpec);
fclose(fid);
CircuitMat{1,4} = double(CircuitMat{1,4});
for RLindex=1:(size(CircuitMat{1,1}))
    if strcmp('RL',CircuitMat{1,1}(RLindex))
        if CircuitMat{1,2}(RLindex) == 0 
            Node1 = 0;
            Node2 = CircuitMat{1,3}(RLindex);
        else
            Node1 = CircuitMat{1,2}(RLindex);
            Node2 = CircuitMat{1,3}(RLindex);
        end
        break
    end
end
Power = [];
RL = [];
if Node1 == 0
    for i=1:1000
    X=modifiednodevoltages(CircuitMat);
    deltaV = X(Node2,1);
    RL(end+1) = CircuitMat{1,4}(RLindex);
    Power(end+1) = ((deltaV).^2)/CircuitMat{1,4}(RLindex);
    CircuitMat{1,4}(RLindex) = CircuitMat{1,4}(RLindex) +0.1;
    end
else
    for i=1:1000
    X=modifiednodevoltages(CircuitMat);
    deltaV = X(Node2,1)-X(Node1,1);
    RL(end+1) = CircuitMat{1,4}(RLindex);
    Power(end+1) = ((deltaV).^2)/CircuitMat{1,4}(RLindex);
    CircuitMat{1,4}(RLindex) = CircuitMat{1,4}(RLindex) +0.1;
    end
end

plot(RL,Power)
[~,powermax] = max(Power);
hold on
plot(RL(powermax),Power(powermax),'or');
title('Power versus Load Resistance');
xlabel('Resistance of RL in ohms');
ylabel('Power dissipated on RL in Watts');
fprintf('Maximum power is %.4fW and corresponding load resistance is %.1f ohms\n',Power(powermax),RL(powermax));

function [x] = modifiednodevoltages(CircuitMat)
%%This functions is same as nodevoltages function but takes a cell array as
%%input instead of a file name and returns an array of node voltages and
%%current passing through voltage source.

%We must group the elements in order to find the number of elements of
%each type
CircuitMat{1,4} = double(CircuitMat{1,4});
[voltagesources,resistances,currentsources] = group(CircuitMat);
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
    if resistances{i,2} ~= 0 && resistances{i,3} ~= 0 
        G(resistances{i,2},resistances{i,3}) = G(resistances{i,2},resistances{i,3})-1/resistances{i,4};
        G(resistances{i,3},resistances{i,2}) = G(resistances{i,3},resistances{i,2})-1/resistances{i,4};
    end
end
%with this for loop we managed to fill the diagonals of G matrix with
%correct values.
for i = 1:n 
    summ = 0;
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
%of the voltage source i is connected to node k then the element B(i,k) =
%1, if negative terminal is connected to k B(i,k) = -1 otherwise it is
%zero.
B = zeros(n,m);
%positive terminals of current and voltage sources are on the third column
%of sourcel cell arrays.
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
