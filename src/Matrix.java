import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.Set;

import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Precision;

public class Matrix {

	private double[][] entries;
	private int numRows, numCols;
	private static final double TOLERANCE = 0.0000000055;   //(1e-8+1e-9)/2..... 1e-9 MAYBE BETTER (1e-11 & "1e-6" common)
	private static final double POLY_TOLERANCE = 0.00001;     //1e-5..... ABSOLUTE BARE MINIMUM (For LaguerreSolver)

	//FUCKING GET RID OF THIS VARIABLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//add in some "final" modifiers (everywhere)

	private int numSwaps = 0;                            //REALLY don't want this to be an instance variable

	public Matrix() {                                    //random square matrix with 1<dim<6 & whole entries -10<x<10
		this.numRows = this.numCols = (int) (4*Math.random()+2);
		this.entries = new double[numRows][numCols];
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < numCols; j++)
				this.entries[i][j] = Math.round(20*Math.random()-10);
		}
	}
	public Matrix(int size) {                            //n*n Identity constructor
		this.numRows = this.numCols = size;
		this.entries = new double[size][size];
		for(int i=0; i < size; i++)
			this.entries[i][i] = 1;
	}
	public Matrix(int numRows, int numCols) {            //m*n Zero matrix constructor
		this.numRows = numRows;
		this.numCols = numCols;
		this.entries = new double[numRows][numCols];
	}
	public Matrix(double[][] entries) {                  //general m*n matrix constructor
		this.numRows = entries.length;
		this.numCols = entries[0].length;
		this.entries = entries;
	}

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
	public void type1eco(int col1, int col2) {
		Matrix result = this.transpose();
		result.type1ero(col1, col2);
		result = result.transpose();
		this.entries = result.entries;
	}
	public void type2eco(double constant, int col) {
		Matrix result = this.transpose();
		result.type2ero(constant, col);
		result = result.transpose();
		this.entries = result.entries;
	}
	public void type3eco(double constant, int col1, int col2) {
		Matrix result = this.transpose();
		result.type3ero(constant, col1, col2);
		result = result.transpose();
		this.entries = result.entries;
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
	public Matrix rightMultiply(Matrix otherMatrix) throws IllegalArgumentException {   //right multiplication

		if(this.numCols != otherMatrix.numRows) 
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		double[][] productEntries = new double[this.numRows][otherMatrix.numCols];
		double[][] otherMatrixEntries = otherMatrix.entries;

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < otherMatrix.numCols; j++) {
				for(int k=0; k < this.numCols; k++)
					productEntries[i][j] += this.entries[i][k]*otherMatrixEntries[k][j];
			}
		}
		Matrix productMatrix = new Matrix(productEntries);
		productMatrix.neutralizeEntries();
		return productMatrix;
	}
	public Matrix scalarMultiply(double scalar) {

		double[][] scaledEntries = this.getEntriesCopy(this.entries);

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				scaledEntries[i][j] *= scalar;
			}
		}
		return new Matrix(scaledEntries);
	}
	public boolean equals(Matrix otherMatrix) throws IllegalArgumentException {   //equals  (AEAP...?)

		if(this.numRows != otherMatrix.numRows && this.numCols != otherMatrix.numCols)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				if(this.entries[i][j] != otherMatrix.entries[i][j]) {    //reverse for efficiency boost...?
					return false;
				}
		}
		return true;
	}

	public int rank() {    //rank (AEAP... better way?)

		int rank = 0;
		Matrix refMatrix = this.toREF();                        //use this to check

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				if(refMatrix.entries[i][j] != 0) {
					rank++;
					break;
				}
			}
		}
		return rank;
	}
	public double trace() throws RuntimeException {

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Inverse not defined.");
		double trace = 0;
		for(int i=0; i < this.numCols; i++)
			trace += this.entries[i][i];
		return trace;
	}
	public double det() throws RuntimeException {     //determinant (AEAP)   +   >>>>>...WRITE A RECURSIVE VERSION...>>>>>

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Determinant not defined.");
		double absDeterminant = 1;
		Matrix refMatrix = this.toREF();    //this initializes numSwaps

		for(int i=0; i < this.numCols; i++) {
			absDeterminant *= refMatrix.entries[i][i];   //using diagonal property of determinants
		}
		double nearestWholeNum = Math.pow(-1,numSwaps)*Math.round(absDeterminant);
		double determinant = Math.pow(-1,numSwaps)*absDeterminant;
		numSwaps = 0;                                                         //KEEP..?

		if(Math.abs(nearestWholeNum-determinant) < Matrix.TOLERANCE)
			return nearestWholeNum;
		return determinant;
	}
	public double[] characteristicPolynomial() {   //(AEAP...?) (Easiest way I know!)

		double[] coefficients = new double[this.numCols+1];     //using secret formuler
		double[] traces = new double[this.numCols+1];
		Matrix originalMatrix = this;
		Matrix tempMatrix = originalMatrix;
		coefficients[this.numCols] = Math.pow(-1,this.numCols);    //Base cases
		traces[this.numCols] = 1;

		for(int i=this.numCols-1; i > -1; i--) {
			traces[i] = tempMatrix.trace();                         //saving space
			tempMatrix = tempMatrix.rightMultiply(originalMatrix);	
			for(int j=this.numCols-1; j > i-1; j--)
				coefficients[i] += (-1.0/(this.numCols-i))*coefficients[this.numCols-j+i]*traces[j];    //reverse formula
			coefficients[i] = Matrix.neutralizeEntry(coefficients[i]);
		}
		return coefficients;
	}
	public LinkedHashMap<Complex,Integer> eigenvalues() {     //complex eigenvalues with respective multiplicites
		
		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Eigenvalues not defined.");

		final double[] coefficients = this.characteristicPolynomial();     //keep finals? apply elsewhere as well..?
		final LaguerreSolver solver = new LaguerreSolver(0.5*TOLERANCE);     //WHY DOESN'T THIS WORK???????
		Complex[] eigenvalues = solver.solveAllComplex(coefficients, 0);
		
		for(int i=0; i < eigenvalues.length; i++) {	
			eigenvalues[i] = new Complex(Matrix.neutralizeValue(eigenvalues[i].getReal()),   //neutralizing
										 Matrix.neutralizeValue(eigenvalues[i].getImaginary()));
		}
		Arrays.sort(eigenvalues, new Comparator<Complex>() {   //sorting the eigenvalues the way I would sort them
			@Override
			public int compare(Complex o1, Complex o2) {
				if(o1.getImaginary() == 0 && o2.getImaginary() != 0) {
					return -1;
				} else if(o1.getImaginary() != 0 && o2.getImaginary() == 0) {
					return 1;
				} else {
					if(o1.getReal()-o2.getReal() < -TOLERANCE)
						return -1;
					else if(Precision.equals(o1.getReal(), o2.getReal(), TOLERANCE))
						if(o1.getImaginary() < o2.getImaginary())
							return -1;
						else
							return 1;
					else
						return 1;
				}
			}
		});
		LinkedHashMap<Complex,Integer> distinctEigenvalues = new LinkedHashMap<Complex,Integer>(this.numCols);
		int multiplicity = 1;
		distinctEigenvalues.put(eigenvalues[0], 1);
		
		for(int i=1; i < this.numCols; i++) {    //takes advantage of the fact that the eigenvalues are sorted
			if(Complex.equals(eigenvalues[i], eigenvalues[i-1], TOLERANCE)) {   //grouping algebraic multiplicities
				eigenvalues[i] = eigenvalues[i-1];      //USE CONTAINS INSTEAD?????????????
				distinctEigenvalues.replace(eigenvalues[i], (++multiplicity)-1, multiplicity);   //better (easier to understand) way...?
			} else {
				multiplicity = 1;
				distinctEigenvalues.put(eigenvalues[i],  multiplicity);
			}
		}
		return distinctEigenvalues;
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

	public Matrix transpose() {   //transpose

		double[][] transposeEntries = new double[this.numCols][this.numRows];

		for(int i=0; i < this.numCols; i++) {
			for(int j=0; j < this.numRows; j++)
				transposeEntries[i][j] = this.entries[j][i]; 
		}
		return new Matrix(transposeEntries);
	}
	public Matrix inverse() throws RuntimeException {    //inverse (AEAP... better way?)

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Inverse not defined.");
		if(this.det() == 0)
			throw new RuntimeException("Matrix has determinant equal to 0. Matrix not invertible.");
		double[][] augmentedEntries = new double[this.numRows][2*this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				augmentedEntries[i][j] = this.entries[i][j];
			augmentedEntries[i][this.numCols+i] = 1;     //saving some time (plugging in identity)
		}
		Matrix rrefAugmentedMatrix = (new Matrix(augmentedEntries)).toRREF();
		double[][] inverseEntries = new double[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				inverseEntries[i][j] = rrefAugmentedMatrix.entries[i][this.numCols+j];   //grabbing 2nd half
		}		
		return new Matrix(inverseEntries);
	}

	/*
	 * //////////////// CONSTRUCTION ZONE BELOW ///////////////////////////////////////////////////////////////////////
	 */

	public ComplexMatrix toD() throws RuntimeException {                      //UNDER CONSTRUCTION

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Diagonal Representation not defined.");

		LinkedHashMap<Complex,Integer> eigenvalues = this.eigenvalues();
		Set<Complex> distinctEigenvalues = eigenvalues.keySet();
		Complex[][] diagonalEntries = new Complex[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				diagonalEntries[i][j] = Complex.ZERO;      //Better way to do this...?
			}
		}
		int count = 0;

		for(Complex eigenvalue: distinctEigenvalues) {
			for(int i=0; i < eigenvalues.get(eigenvalue); i++)
				diagonalEntries[count][count++] = eigenvalue;
		}		
		return new ComplexMatrix(diagonalEntries);
	}


	/*
	 * TODO: MAKE reduceToREF & reduceToRREF (void) and/or Singular Value Decomposition (SVD)
	 */


	/*
	 * ///////////// "BEHIND-THE-SCENES METHODS" BELOW ////////////////////////////////////////////////////////////////
	 */

	public static double neutralizeValue(double entry) {
		double nearestWholeNum = Math.round(entry);
		if(Math.abs(nearestWholeNum-entry) < Matrix.POLY_TOLERANCE)    //MAKE THIS DISTINCT (NAMING)
			return nearestWholeNum;
		return entry;
	}
	private static double neutralizeEntry(double entry) {
		double nearestWholeNum = Math.round(entry);
		if(Math.abs(nearestWholeNum-entry) < Matrix.TOLERANCE)     //MAKE THIS DISTINCT (NAMING)
			return nearestWholeNum;
		return entry;
	}
	private void neutralizeEntries() {
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				/*
				double nearestWholeNum = Math.round(this.entries[i][j]);
				if(Math.abs(nearestWholeNum-this.entries[i][j]) < Matrix.TOLERANCE)
					this.entries[i][j] = nearestWholeNum;
				 */
				this.entries[i][j] = Matrix.neutralizeEntry(this.entries[i][j]);
			}
		}
	}

	//TODO: MAKE SOME SETTERS...?

	public double[][] getEntriesCopy(double[][] entries) {
		double[][] entriesCopy = new double[entries.length][entries[0].length];
		for(int i=0; i < entries.length; i++) {
			for(int j=0; j < entries[0].length; j++)     //ACTUALLLY saving the original entries
				entriesCopy[i][j] = entries[i][j];
		}
		return entriesCopy;
	}

	public double[][] getEntries() {                 //get matrix content (array)   
		return this.entries; 
	}
	public int getNumRows() {                              //get number of rows
		return this.numRows;
	}
	public int getNumCols() {                              //get number of columns
		return this.numCols;
	}
	
	/*
	 * ///////////// "DEBUGGING METHODS" BELOW ////////////////////////////////////////////////////////////////////////
	 */

	public static double evaluate(double[] coefficients, double value) {
		double result = 0;
		for(int i=0; i < coefficients.length; i++) {
			result += coefficients[i]*Math.pow(value,i);
		}
		return result;
	}

	public static void printComplexArray(Complex[] roots) {
		String result = "[";
		for(int i=0; i < roots.length; i++) {
			result += roots[i].getReal() + " + " + roots[i].getImaginary() + "i, ";
		}
		result += "]";
		System.out.println(result);
	}
	public static void printArray(double[] entries) {
		String result = "[";
		for(int i=0; i < entries.length; i++) {
			result += entries[i] + ", ";
		}
		result += "]";
		System.out.println(result);
	}
	public static void print2DArray(double[][] entries) {
		String result = "[";
		for(int i=0; i < entries.length; i++) {
			result += "[";
			for(int j=0; j < entries[0].length; j++)
				result += entries[i][j] + ", ";
			result += "]";
		}
		result += "]";
		System.out.println(result);
	}

	@Override
	public String toString() {                       //basic toString
		this.neutralizeEntries();
		String result = "";
		for(int i=0; i < this.numRows; i++) {
			result += "| ";
			for(int j=0; j < this.numCols; j++)
				result += String.format("% 8.3f", this.entries[i][j]) + " ";
			result += " |\n";
		}
		return result;
	}

	public String toStringLong() {              //special toString
		this.neutralizeEntries();
		String result = "";
		for(int i=0; i < this.numRows; i++) {
			result += "| ";
			for(int j=0; j < this.numCols; j++)
				result += String.format("% 9.2f", this.entries[i][j]) + " ";
			result += " |\n";
		}
		return result;
	}

	/*
	 * TODO: MAKE TOSTRING(s): SMALL NUMS (floats < 1) & BIG NUMS (longs > 1000)
	 */

}