
# Simple Circuit Simulator Using Modified Node Analysis (MNA)

A simple circuit simulator that is capable of simulating circuits with the circuit elements:

- Resistors
- Independent Voltage Sources
- Independent Current Sources

Program reads an input txt file representing the circuit and determines the node voltages.

The text file will include the information regarding circuits with the following rules:
- Each element is entered in a single row.
- The first column is the unique identifier for the element whose first letter indicates the type of the
element: R, I or V and the rest is an integer. The second and the third columns denote the node numbers
of the element. The last column denotes the value of the element in Ohms, Amperes or Volts.
- Node Number at the Second Column < Node Number at the Third Column.
- Positive value for the current source means that the current is entering the node at the Third Column.
- Positive value for the voltage source means: Voltage of node at the Second Column < Voltage of node at the Third Column.
- For example the following circuit will be entered as follows:

![image](https://user-images.githubusercontent.com/63296692/153558428-20ac80a0-2679-4a97-aea2-d86afce091af.png)
![image](https://user-images.githubusercontent.com/63296692/153558387-8a02a3bf-a7f3-46ed-96d1-f5e4d55c157c.png)


### MNA Algorithm [^1]

MNA applied to a circuit with only passive elements (resistors) and independent current and voltage sources resultsin a matrix equation of the form:

A ∙ x = z


For a circuit with n nodes and m independent voltage sources:
- The A matrix, stated above:
  - is (n+m)x(n+m) in size, and consists only of known quantities.
  - the nxn part of the matrix in the upper left:
     - has only passive elements
     - elements connected to ground appear only on the diagonal
     - elements not connected to ground are both on the diagonal and off-diagonal terms.
- The rest of the A matrix (not included in the nxn upper left part) contains only 1, -1 and 0.
- The x matrix:
  - is an (n+m)x1 vector that holds the unknown quantities (node voltages and the currents through the independent voltage sources).
  - the top n elements are the n node voltages.
  - the bottom m elements represent the currents through the m independent voltage sources in the circuit.
- The z matrix:
  - is an (n+m)x1 vector that holds only known quantities
  - the top n elements are either zero or the sum and difference of independent current sources in the circuit.
  - the bottom m elements represent the m independent voltage sources in the circuit.
The unknown quantities can be obtained by solving the linear system of equations A ∙ x = z .

### Power_wrt_load

Another program that uses the program above as a function in a loop to sweep the load resistance, then calculates and plots the power (W) dissipated in the load versus the load resistance. Also displays the maximum power and the corresponding resistance on command window as well as puts a marker on the plot.


[^1]: http://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA3.html
