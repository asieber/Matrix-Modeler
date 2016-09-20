package model;
import java.util.LinkedHashSet;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Precision;

/*
 * The ComplexMatrix class models matrices with complex entries and supports basic matrix operations.
 */

public class ComplexMatrix {

	private Complex[][] entries;
	private int numRows, numCols;
	private static final double ROUND_TOLERANCE = 0.0000000001;
	
	private int numSwaps;

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
	public ComplexMatrix(double[][] realEntries) {          //CHANGE ABOVE TO "realEntries"....!!!!!!!!!!!!!!!!!!!
		this.numRows = realEntries.length;
		this.numCols = realEntries[0].length;
		this.entries = new Complex[this.numRows][this.numCols];
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				this.entries[i][j] = new Complex(realEntries[i][j], 0);
		}
	}
	
	public void type1ero(int row1, int row2) {
		Complex[] tempRowEntries = this.entries[row1-1];
		this.entries[row1-1] = this.entries[row2-1];
		this.entries[row2-1] = tempRowEntries;
	}
	public void type2ero(Complex constant, int row) {
		for(int i=0; i < this.numCols; i++)
			this.entries[row-1][i] = constant.multiply(this.entries[row-1][i]); 
		this.roundEntries();
	}
	public void type3ero(Complex constant, int row1, int row2) {
		for(int i=0; i < this.numCols; i++)
			this.entries[row2-1][i] = this.entries[row2-1][i].add(constant.multiply(this.entries[row1-1][i])); 
		this.roundEntries();
	}
	
	//MAKE A METHOD TO MULTIPLY BY A REAL MATRIX (& OTHERS INVOLVING REAL MATRICES/VALUES)

	public ComplexMatrix plus(ComplexMatrix otherMatrix) throws IllegalArgumentException {   //plus

		if(this.numRows != otherMatrix.numRows || this.numCols != otherMatrix.numCols)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		Complex[][] sumEntries = new Complex[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				sumEntries[i][j] = this.entries[i][j].add(otherMatrix.entries[i][j]);
		}
		return new ComplexMatrix(sumEntries);
	}
	public ComplexMatrix minus(ComplexMatrix otherMatrix) throws IllegalArgumentException {   //minus

		if(this.numRows != otherMatrix.numRows || this.numCols != otherMatrix.numCols)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		Complex[][] sumEntries = new Complex[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				sumEntries[i][j] = this.entries[i][j].subtract(otherMatrix.entries[i][j]);
		}
		return new ComplexMatrix(sumEntries);
	}
	public ComplexMatrix rightMultiply(ComplexMatrix otherMatrix) throws IllegalArgumentException {   //right multiplication

		if(this.numCols != otherMatrix.numRows) 
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		Complex[][] productEntries = new Complex[this.numRows][otherMatrix.numCols];
		
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				productEntries[i][j] = Complex.ZERO;    //EXTRA LOOP NECESSARY FOR COMPLEX (but not real).....?
		}
		Complex[][] otherMatrixEntries = otherMatrix.entries;

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < otherMatrix.numCols; j++) {
				for(int k=0; k < this.numCols; k++)
					productEntries[i][j] = productEntries[i][j].add(this.entries[i][k].multiply(otherMatrixEntries[k][j]));
			}
		}
		ComplexMatrix productMatrix = new ComplexMatrix(productEntries);
		productMatrix.roundEntries();
		return productMatrix;
	}
	public ComplexMatrix scalarMultiply(double scalar) {

		Complex[][] scaledEntries = this.getEntriesCopy(this.entries);   //MAKE A COMPLEX COPIER

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				scaledEntries[i][j] = scaledEntries[i][j].multiply(scalar);
			}
		}
		return new ComplexMatrix(scaledEntries);
	}
	public ComplexMatrix scalarMultiply(Complex scalar) {

		Complex[][] scaledEntries = this.getEntriesCopy(this.entries);   //MAKE A COMPLEX COPIER

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				scaledEntries[i][j] = scaledEntries[i][j].multiply(scalar);
			}
		}
		return new ComplexMatrix(scaledEntries);
	}
	public ComplexMatrix toPower(int power) {
		
		ComplexMatrix originalMatrix = this;
		ComplexMatrix result = originalMatrix;
		
		for(int i=0; i < power-1; i++) {
			result = result.rightMultiply(originalMatrix);
		}
		return result;
	}
	public boolean equals(ComplexMatrix otherMatrix) throws IllegalArgumentException {   //equals  (AEAP...?)

		if(this.numRows != otherMatrix.numRows && this.numCols != otherMatrix.numCols)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				if(!this.entries[i][j].equals(otherMatrix.entries[i][j])) {    //reverse for efficiency boost...?
					return false;
				}
		}
		return true;
	}
	
	public int rank() {    //rank (AEAP... better way?)

		int rank = 0;
		ComplexMatrix refMatrix = this.toREF();                        //use this to check

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				if(!refMatrix.entries[i][j].equals(Complex.ZERO)) {
					rank++;
					break;
				}
			}
		}
		return rank;
	}
	public Complex det() throws RuntimeException {     //determinant (AEAP)   +   >>>>>...WRITE A RECURSIVE VERSION...>>>>>

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Determinant not defined.");
		Complex absDeterminant = Complex.ONE;
		ComplexMatrix refMatrix = this.toREF();    //this initializes numSwaps

		for(int i=0; i < this.numCols; i++) {
			absDeterminant = absDeterminant.multiply(refMatrix.entries[i][i]);   //using diagonal property of determinants
		}
		
		//GOTTA WORK ON "ROUND" FOR COMPLEX
		
		//Complex nearestWholeNum = Math.pow(-1,numSwaps)*Math.round(absDeterminant);
		
		Complex determinant = absDeterminant.multiply(Math.pow(-1,numSwaps));
		numSwaps = 0;                                                         //KEEP..?
		
		//GOTTA WORK ON "ROUND" FOR COMPLEX

		//if(Math.abs(nearestWholeNum-determinant) < Matrix.TOLERANCE)
		//	return nearestWholeNum;
		
		return determinant;
	}

	public ComplexMatrix toREF () {   //row echelon form

		Complex[][] originalEntries = this.getEntriesCopy(this.entries);
		int startRowIndex = 0;
		numSwaps = 0;

		for(int i=0; i < this.numCols; i++) {
			if(clearColumnDown(i, startRowIndex))   //if the column is nonzero, pivot advances to next row
				startRowIndex++; 
		}
		ComplexMatrix refMatrix = new ComplexMatrix(this.entries);
		this.entries = originalEntries;                  //setting matrix back to its original form)
		return refMatrix;
	}
	private boolean clearColumnDown(int colIndex, int startRowIndex) {    //helper method (toREF)

		boolean isNonZeroCol = false;

		for(int i=startRowIndex; i < this.numRows; i++) {
			if(!this.entries[i][colIndex].equals(Complex.ZERO)) {
				if(i != startRowIndex)
					numSwaps++;                            //BETTER WAY TO DO THIS?????????????
				this.type1ero(i+1, startRowIndex+1); 
				isNonZeroCol = true;
				break;
			}
		}
		if(isNonZeroCol) {
			for(int i = startRowIndex; i < this.numRows-1; i++)
				this.type3ero(this.entries[i+1][colIndex].multiply(-1.0).divide(this.entries[startRowIndex][colIndex]), 
								startRowIndex+1, i+2);
		}
		return isNonZeroCol;
	}
	public ComplexMatrix toRREF() {    //reduced row echelon form

		Complex[][] originalEntries = this.getEntriesCopy(this.entries);
		ComplexMatrix refMatrix = this.toREF();
		Complex[][] refEntries = refMatrix.entries;
		int startRowIndex = this.numRows-1;

		for(int i=startRowIndex; i > -1; i--) {
			for(int j=0; j < this.numCols; j++) {
				if(!refEntries[i][j].equals(Complex.ZERO)) {
					refMatrix.type2ero(Complex.ONE.divide(refEntries[i][j]), i+1);
					refMatrix.clearColumnUp(j, i);                    //ordering of args...????
					break;
				}
			}
		}
		ComplexMatrix rrefMatrix = refMatrix;
		this.entries = originalEntries;                  //setting matrix back to its original form
		return rrefMatrix;
	}
	private void clearColumnUp(int colIndex, int startRowIndex) {    //helper method (toRREF)

		for(int i=startRowIndex-1; i > -1; i--) {
			Complex colEntry = this.entries[i][colIndex];
			if(!colEntry.equals(Complex.ZERO))
				this.type3ero(colEntry.multiply(-1.0), startRowIndex+1, i+1);   //no scaling needed, pivot already == 1 from toRREF
		}
	}
	
	public ComplexMatrix transpose() {   //transpose

		Complex[][] transposeEntries = new Complex[this.numCols][this.numRows];

		for(int i=0; i < this.numCols; i++) {
			for(int j=0; j < this.numRows; j++)
				transposeEntries[i][j] = this.entries[j][i]; 
		}
		return new ComplexMatrix(transposeEntries);
	}
	public ComplexMatrix inverse() throws RuntimeException {    //inverse (AEAP... better way?)

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Inverse not defined.");
		if(this.det().equals(Complex.ZERO))
			throw new RuntimeException("Matrix has determinant equal to 0. Matrix not invertible.");
		Complex[][] augmentedEntries = new Complex[this.numRows][2*this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				augmentedEntries[i][j] = this.entries[i][j];
			//augmentedEntries[i][this.numCols+i] = Complex.ONE;     //saving some time (plugging in identity)
		}
		for(int i=0; i < this.numRows; i++) {
			for(int j=this.numCols; j < 2*this.numCols; j++) {
				if(j == i+this.numCols)
					augmentedEntries[i][j] = Complex.ONE;    //EXTRA LOOPS NECESSARY FOR COMPLEX (but not real).....?
				else
					augmentedEntries[i][j] = Complex.ZERO;
			}
		}
		ComplexMatrix rrefAugmentedMatrix = (new ComplexMatrix(augmentedEntries)).toRREF();
		Complex[][] inverseEntries = new Complex[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				inverseEntries[i][j] = rrefAugmentedMatrix.entries[i][this.numCols+j];   //grabbing 2nd half
		}		
		return new ComplexMatrix(inverseEntries);
	}
	
	public static LinkedHashSet<Complex[]> eigenvectors(ComplexMatrix matrix, Complex eigenvalue) {
		
		ComplexMatrix diffMatrix = matrix.minus(new ComplexMatrix(matrix.numCols).scalarMultiply(eigenvalue));
		diffMatrix = diffMatrix.toRREF();
		
		//use rank as shortcut...?
		
		//do some algebra
		
		return null;
		
	}

	public static Complex roundEntry(Complex entry) {
		double originalReal = entry.getReal();
		double originalImag = entry.getImaginary();
		double nearestWholeNumReal = Math.round(entry.getReal());
		double nearestWholeNumImag = Math.round(entry.getImaginary());
		
		//ADD IN: IF NEITHER IS TRUE, RETURN ORIGINAL ENTRY
		
		//JUST USE PRE-EXISTING Complex.equals(..., ...) &&/|| Precision.equals(..., ...)...

		if(Math.abs(nearestWholeNumReal-originalReal) < ComplexMatrix.ROUND_TOLERANCE)     //MAKE THIS DISTINCT (NAMING)
			originalReal = nearestWholeNumReal;
		if(Math.abs(nearestWholeNumImag-originalImag) < ComplexMatrix.ROUND_TOLERANCE)     //MAKE THIS DISTINCT (NAMING)
			originalImag = nearestWholeNumImag;
		return new Complex(originalReal, originalImag);
	}
	public void roundEntries() {
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				
				//double nearestWholeNum = Math.round(this.entries[i][j]);
				//if(Math.abs(nearestWholeNum-this.entries[i][j]) < Matrix.TOLERANCE)
				//	this.entries[i][j] = nearestWholeNum;
				
				this.entries[i][j] = ComplexMatrix.roundEntry(this.entries[i][j]);
			}
		}
	}

	public Complex[][] getEntriesCopy(Complex[][] entries) {
		Complex[][] entriesCopy = new Complex[entries.length][entries[0].length];
		for(int i=0; i < entries.length; i++) {
			for(int j=0; j < entries[0].length; j++)     //ACTUALLLY saving the original entries
				entriesCopy[i][j] = entries[i][j];
		}
		return entriesCopy;
	}
	public Complex[][] getEntriesCopy() {
		Complex[][] entriesCopy = new Complex[this.entries.length][this.entries[0].length];
		for(int i=0; i < this.entries.length; i++) {
			for(int j=0; j < this.entries[0].length; j++)     //ACTUALLLY saving the original entries
				entriesCopy[i][j] = this.entries[i][j];
		}
		return entriesCopy;
	}
	
	public Complex[][] getEntries() {                 //get matrix content (array)   
		return this.entries; 
	}
	public int getNumRows() {                              //get number of rows
		return this.numRows;
	}
	public int getNumCols() {                              //get number of columns
		return this.numCols;
	}

//	DEBUG toString W/:
//
//	A toString: 
//	|    6.000    1.000   -6.000    4.000  |
//	|   -5.000   -4.000    0.000   -1.000  |
//	|    5.000    9.000    6.000   -3.000  |
//	|    0.000   -7.000    6.000    1.000  |
	
	@Override
	public String toString() {                       //basic toString
		
		this.roundEntries();

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
