import java.util.Hashtable;
import java.util.Set;

import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Precision;

public class ComplexMatrix {


	private Complex[][] entries;
	private int numRows, numCols;
	private static final double TOLERANCE = 0.0000000055;

	public ComplexMatrix() {                        //random square ComplexMatrix with 1<dim<6 & whole entries -10<x<10
		this.numRows = this.numCols = (int) (4*Math.random()+2);
		this.entries = new Complex[numRows][numCols];
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < numCols; j++)
				this.entries[i][j] = new Complex(Math.round(20*Math.random()-10), 
						Math.round(20*Math.random()-10));
		}
	}
	public ComplexMatrix(int size) {                            //n*n Identity constructor
		this.numRows = this.numCols = size;
		this.entries = new Complex[size][size];
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				if(i == j)
					this.entries[i][j] = Complex.ONE;
				else
					this.entries[i][j] = Complex.ZERO;
			}
		}
	}
	public ComplexMatrix(int numRows, int numCols) {            //m*n Zero ComplexMatrix constructor
		this.numRows = numRows;
		this.numCols = numCols;
		this.entries = new Complex[numRows][numCols];
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				this.entries[i][j] = Complex.ZERO;
			}
		}
	}
	public ComplexMatrix(Complex[][] complexEntries) {          //CHANGE ABOVE TO "realEntries"....!!!!!!!!!!!!!!!!!!!
		this.numRows = complexEntries.length;
		this.numCols = complexEntries[0].length;
		this.entries = complexEntries;
	}

	/*
	public void type1ero(int row1, int row2) {
		double[] tempRowEntries = this.entries[row1-1];
		this.entries[row1-1] = this.entries[row2-1];
		this.entries[row2-1] = tempRowEntries;
	}
	public void type2ero(double constant, int row) {
		for(int i=0; i < this.numCols; i++)
			this.entries[row-1][i] *= constant; 
		this.neutralizeEntries();
	}
	public void type3ero(double constant, int row1, int row2) {
		for(int i=0; i < this.numCols; i++)
			this.entries[row2-1][i] += constant*this.entries[row1-1][i]; 
		this.neutralizeEntries();
	}

	public Matrix plus(Matrix otherMatrix) throws IllegalArgumentException {   //plus

		if(this.numRows != otherMatrix.numRows || this.numCols != otherMatrix.numCols)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		double[][] sumEntries = new double[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				sumEntries[i][j] = this.entries[i][j]+otherMatrix.entries[i][j];
		}
		return new Matrix(sumEntries);
	}
	public Matrix minus(Matrix otherMatrix) throws IllegalArgumentException {   //minus

		if(this.numRows != otherMatrix.numRows || this.numCols != otherMatrix.numCols)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		double[][] sumEntries = new double[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				sumEntries[i][j] = this.entries[i][j]-otherMatrix.entries[i][j];
		}
		return new Matrix(sumEntries);
	}
	public Matrix scalarMultiply(Complex scalar) {

		Complex[][] scaledEntries = null; //this.getEntriesCopy(this.entries);   //MAKE A COMPLEX COPIER

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				scaledEntries[i][j] = scaledEntries[i][j].multiply(scalar);
			}
		}
		return new Matrix(scaledEntries);
	}

	public Matrix toREF () {   //row echelon form

		double[][] originalEntries = this.getEntriesCopy(this.entries);
		int startRowIndex = 0;
		numSwaps = 0;

		for(int i=0; i < this.numCols; i++) {
			if(clearColumnDown(i, startRowIndex))   //if the column is nonzero, pivot advances to next row
				startRowIndex++; 
		}
		Matrix refMatrix = new Matrix(this.entries);
		this.entries = originalEntries;                  //setting matrix back to its original form)
		return refMatrix;
	}
	private boolean clearColumnDown(int colIndex, int startRowIndex) {    //helper method (toREF)

		boolean isNonZeroCol = false;

		for(int i=startRowIndex; i < this.numRows; i++) {
			if(this.entries[i][colIndex] != 0) {
				if(i != startRowIndex)
					numSwaps++;                            //BETTER WAY TO DO THIS?????????????
				this.type1ero(i+1, startRowIndex+1); 
				isNonZeroCol = true;
				break;
			}
		}
		if(isNonZeroCol) {
			for(int i = startRowIndex; i < this.numRows-1; i++)
				this.type3ero(-this.entries[i+1][colIndex]/this.entries[startRowIndex][colIndex], startRowIndex+1, i+2);
		}
		return isNonZeroCol;
	}
	public Matrix toRREF() {    //reduced row echelon form

		double[][] originalEntries = this.getEntriesCopy(this.entries);
		Matrix refMatrix = this.toREF();
		double[][] refEntries = refMatrix.entries;
		int startRowIndex = this.numRows-1;

		for(int i=startRowIndex; i > -1; i--) {
			for(int j=0; j < this.numCols; j++) {
				if(refEntries[i][j] != 0) {
					refMatrix.type2ero(1/refEntries[i][j], i+1);
					refMatrix.clearColumnUp(j, i);                    //ordering of args...????
					break;
				}
			}
		}
		Matrix rrefMatrix = refMatrix;
		this.entries = originalEntries;                  //setting matrix back to its original form
		return rrefMatrix;
	}
	private void clearColumnUp(int colIndex, int startRowIndex) {    //helper method (toRREF)

		for(int i=startRowIndex-1; i > -1; i--) {
			double colEntry = this.entries[i][colIndex];
			if(colEntry != 0)
				this.type3ero(-colEntry, startRowIndex+1, i+1);   //no scaling needed, pivot already == 1 from toRREF
		}
	}
	 */

	/*
	public Matrix diagonalize() {                      //UNDER CONSTRUCTION

		final double[] coefficients = this.characteristicPolynomial();     //keep finals? apply elsewhere as well..?
		final LaguerreSolver solver = new LaguerreSolver(0.5*TOLERANCE);     //WHY DOESN'T THIS WORK???????
		Complex[] eigenvalues = solver.solveAllComplex(coefficients, 0);

		for(int i=0; i < eigenvalues.length; i++) {	
			eigenvalues[i] = new Complex(Matrix.neutralizeValue(eigenvalues[i].getReal()),   //neutralizing
										 Matrix.neutralizeValue(eigenvalues[i].getImaginary()));
		}
		for(int i=1; i < eigenvalues.length; i++) {     //making similar eigenvalue(s) identical (NECESSARY?)
			for(int j=i-1; j > -1; j--) {
				if(Complex.equals(eigenvalues[i], eigenvalues[i-1], TOLERANCE))
					eigenvalues[i] = eigenvalues[j];
				if(Precision.equals(eigenvalues[i].getReal(), -eigenvalues[j].getReal(), TOLERANCE))
					eigenvalues[i] = new Complex(-eigenvalues[j].getReal(), eigenvalues[i].getImaginary());
				if(Precision.equals(eigenvalues[i].getReal(), eigenvalues[j].getReal(), TOLERANCE))
					eigenvalues[i] = new Complex(eigenvalues[j].getReal(), eigenvalues[i].getImaginary());
				if(Precision.equals(eigenvalues[i].getImaginary(), -eigenvalues[j].getImaginary(), TOLERANCE))
					eigenvalues[i] = new Complex(eigenvalues[i].getReal(), -eigenvalues[j].getImaginary());
				if(Precision.equals(eigenvalues[i].getImaginary(), eigenvalues[j].getImaginary(), TOLERANCE))
					eigenvalues[i] = new Complex(eigenvalues[i].getReal(), eigenvalues[j].getImaginary());
			}
		}

		System.out.println("Eigenvalues: ");
		Matrix.printComplexArray(eigenvalues);      //DEBUGGING
		System.out.println();

		Hashtable<Complex,Integer> table = new Hashtable<Complex,Integer>(this.numCols);
		int multiplicity = 1;
		table.put(eigenvalues[0], 1);

		for(int i=1; i < this.numCols; i++) {
			if(Complex.equals(eigenvalues[i], eigenvalues[i-1], TOLERANCE)) {   //grouping algebraic multiplicities
				eigenvalues[i] = eigenvalues[i-1];
				//multiplicity++;
				table.replace(eigenvalues[i], (++multiplicity)-1, multiplicity);   //better (easier to understand) way...?
			} else {
				multiplicity = 1;
				table.put(eigenvalues[i],  multiplicity);
			}
		}
		Set<Complex> distinctEigenvalues = table.keySet();
		Matrix eigenMatrix = null;

		for(Complex eigenvalue: distinctEigenvalues) {

			System.out.println("Multiplicity of (" + eigenvalue.getReal() + " + " + eigenvalue.getImaginary() 
								+ "i) is: " + table.get(eigenvalue));


			//eigenMatrix = this.minus(new Matrix(this.numCols).scalarMultiply(eigenvalue));
			//Matrix rrefEigenMatrix = eigenMatrix.toRREF();     //gonna need to define this for complex entries


			//DO SOME ALGEBRA TO FIGURE OUT EIGENVECTORS.....
		}
		//error checking if diagonalizable or not (use method: Matrix.eigenspaceDimension(Matrix otherMatrix))

		return null;
	}
	 */

	@Override
	public String toString() {                       //basic toString
		
		//this.neutralizeEntries();

		String result, real, imagNum, imag;
		result = real = imagNum = imag = "";

		for(int i=0; i < this.numRows; i++) {
			result += "| ";
			for(int j=0; j < this.numCols; j++) {
				real = String.format("% 8.3f", this.entries[i][j].getReal());
				imagNum = String.format("%.3f", Math.abs(this.entries[i][j].getImaginary())) + "i";
				imag = String.format("%-8s", imagNum);				
				result += real;
				if(this.entries[i][j].getImaginary() >= 0)
					result += " + ";
				else
					result += " - ";
				result += imag;
			}
			result += " |\n";
		}
		return result;
	}

}
