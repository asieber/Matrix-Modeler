package model;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Precision;

/*
 * The Matrix class models matrices with real entries and supports basic matrix operations.
 */

public class Matrix {

	private double[][] entries;
	private int numRows, numCols;
	
	//MAKE SURE POLY_TOLERANCE & ROUND_TOLERANCE ARE BEING USED PROPERLY AND LOGICALLY THROUGHOUT CODE.....
	
	/*ROUND_TOLERANCE is used to round individual terms to their true nearest whole number values. It is also used to
	allow entries to be significantly small without becoming zero as a result of inaccuracies caused by
	algebraic manipulation. ROUND_TOLERANCE is a general tolerance for matrix entries and is used for individual 
	values whereas POLY_TOLERANCE is typically pairwise and is strictly used in special cases pertaining to the 
	LaguerreSolver.*/
	
	private static final double ROUND_TOLERANCE = 0.0000000001; /*.01*POLY_TOLERANCE;*/   //SWITCH TO (1e-8+1e-9)/2...?
	
	/*POLY_TOLERANCE is used to equate similar terms that are truly equal but came out of the LaguerreSolver
	unequal. It is also used to round terms to their true nearest whole number values if they
	came out of the LaguerreSolver with an error of about +/- 0.001 or less. ROUND_TOLERANCE
	is used for the same purpose as the latter however it is not always small enough to be used 
	for the values obtained from the LaguerreSolver. Problems really only arose after trying to obtain the 
	eigenvalues of whole-entry diagonal matrices. Certain eigenvalues that should have been equal appeared to have been 
	coming out within about +/- .001 of each other and were also about 0.001 away from their true nearest whole
	number value (which is a pretty significant difference).*/
	
	private static final double LAGUERRE_TOLERANCE = 0.001; /*0.00001;*/   //(For LaguerreSolver & MANY OTHERS NOW: eigenvalues, roundComplexEntry, etc...)
	
	//add in some "final" modifiers (everywhere)

	//GET RID OF THIS VARIABLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	private int numSwaps = 0;                            //REALLY don't want this to be an instance variable
	
	public Matrix() {                                    //random square matrix with 1<dim<6 & whole entries -10<x<10
		this.numRows = this.numCols = /*2;*/ (int) (4*Math.random()+2);
		this.entries = new double[numRows][numCols];
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < numCols; j++)
				this.entries[i][j] = Math.round(20*Math.random()-10);
		}
		//this.neutralizeEntries();
	}
	public Matrix(Random rand) {                         //dep injected matrix with 1<dim<6 & whole entries -10<x<10
		this.numRows = this.numCols = /*2;*/ (int) (4*rand.nextFloat()+2);
		this.entries = new double[numRows][numCols];
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < numCols; j++)
				this.entries[i][j] = Math.round(20*rand.nextFloat()-10);
		}
		//this.neutralizeEntries();
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
		//this.neutralizeEntries();
	}

	public void type1ero(int row1, int row2) {
		double[] tempRowEntries = this.entries[row1-1];
		this.entries[row1-1] = this.entries[row2-1];
		this.entries[row2-1] = tempRowEntries;
	}
	public void type2ero(double constant, int row) {
		for(int i=0; i < this.numCols; i++)
			this.entries[row-1][i] *= constant; 
		this.roundEntries();
	}
	public void type3ero(double constant, int row1, int row2) {
		for(int i=0; i < this.numCols; i++)
			this.entries[row2-1][i] += constant*this.entries[row1-1][i]; 
		this.roundEntries();
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

	//MAKE A METHOD TO MULTIPLY BY A COMPLEX MATRIX (THINK OF OTHER INVOLVING COMPLEX VALUES (plus, minus, etc...!!!)
	
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
		productMatrix.roundEntries();
		return productMatrix;
	}
	public Matrix scalarMultiply(double scalar) {

		double[][] scaledEntries = this.getEntriesCopy();

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++) {
				scaledEntries[i][j] *= scalar;
			}
		}
		return new Matrix(scaledEntries);
	}
	public ComplexMatrix scalarMultiply(Complex scalar) {
		
		double[][] originalEntries = this.getEntriesCopy();
		Complex[][] scaledEntries = new Complex[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				scaledEntries[i][j] = scalar.multiply(originalEntries[i][j]);
		}
		return new ComplexMatrix(scaledEntries);
	}
	public Matrix toPower(int power) {
		
		Matrix originalMatrix = this;
		Matrix result = originalMatrix;
		
		for(int i=0; i < power-1; i++) {
			result = result.rightMultiply(originalMatrix);
		}
		return result;
	}
	public boolean equals(Matrix otherMatrix) throws IllegalArgumentException {   //equals  (AEAP...?)

		if(this.numRows != otherMatrix.numRows && this.numCols != otherMatrix.numCols)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				if(this.entries[i][j] != otherMatrix.entries[i][j])    //reverse for efficiency boost...?
					return false;
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
			throw new RuntimeException("Matrix isn't square. Trace not defined.");
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

		if(Math.abs(nearestWholeNum-determinant) < Matrix.ROUND_TOLERANCE)
			return nearestWholeNum;
		return determinant;
	}
	public double[] characteristicPolynomial() throws RuntimeException {  //(AEAP...?) (Easiest way I know!)

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Characteristic polynomial not defined.");
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
			coefficients[i] = Matrix.roundEntry(coefficients[i]);    //neutralizeValue...?
		}
		return coefficients;
	}
	
	public LinkedHashSet<double[]> range() {                 //JUST USE Set INSTEAD OF LinkedHashSet????????????????
		
		Matrix rrefMatrix = this.toRREF();              //USE copy "METHOD"...???????????
		double[][] rrefEntries = rrefMatrix.entries;
		Set<Integer> linIndColIndices = new LinkedHashSet<Integer>(this.rank());
		
		for(int i=0; i < rrefMatrix.numRows; i++) {
			for(int j=0; j < rrefMatrix.numCols; j++) {
				if(rrefEntries[i][j] != 0) {
					linIndColIndices.add(j);
					break;
				}
			}
		}
		Matrix transpose = this.transpose();
		double[][] transposeEntries = transpose.entries;                                //JUST USE .entries??????
		LinkedHashSet<double[]> rangeBasis = new LinkedHashSet<double[]>(this.rank());
		
		for(int i=0; i < transpose.numRows; i++) {
			if(linIndColIndices.contains(i)) {
				rangeBasis.add(transposeEntries[i]);
			}
		}
		return rangeBasis;
	}
	public LinkedHashSet<double[]> nullspace() {
		
		double[][] nullspaceEntries = new double[this.numCols][this.numCols];
		Matrix rrefMatrix = new Matrix(this.getEntriesCopy()).toRREF();      //EASIER WAY...?
		double[][] rrefEntries = rrefMatrix.entries;
		
		for(int i=0; i < Math.min(rrefEntries.length, nullspaceEntries.length); i++) {
			nullspaceEntries[i] = rrefEntries[i];
		}
//		for(int i=0; i < nullspaceEntries.length-rrefEntries.length; i++) {    //for Complex/ev only, already 0 for doubles
//			nullspaceEntries[i] = new double[this.numCols];
//		}

		for(int i=0; i < nullspaceEntries.length; i++) {
			for(int j=0; j < nullspaceEntries[i].length; j++) {
				if(j==i) {
					if(nullspaceEntries[j][j] == 0) {
						nullspaceEntries[j][j] = 1;
					} else {
						nullspaceEntries[j][j] = 0;
					}
				} else {
					nullspaceEntries[i][j] *= -1;
				}
			}
		}
		Matrix nullSpaceMatrix = new Matrix(nullspaceEntries);    //nullSpace
		nullSpaceMatrix = nullSpaceMatrix.transpose();
		nullSpaceMatrix.roundEntries();                                  //EVEN NECESSARY...?
		double[][] eigenSpace = nullSpaceMatrix.getEntriesCopy();      //JUST USE entries INSTEAD????????????????????
		LinkedHashSet<double[]> eigenvectors = new LinkedHashSet<double[]>(nullSpaceMatrix.rank());
		
		for(int i=0; i < eigenSpace.length; i++) {
			for(int j=0; j < eigenSpace[i].length; j++) {
				if(eigenSpace[i][j] != 0) {
					eigenvectors.add(eigenSpace[i]);
					break;
				}
			}
		}
		return eigenvectors;    //DOES THIS WORK FOR NON-SQUARE MATRICES...?
	}
	public LinkedHashMap<Complex,Integer> eigenvalues() throws RuntimeException {     //complex eigenvalues with respective multiplicites
		
		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Eigenvalues not defined.");
		
		final double[] coefficients = this.characteristicPolynomial();     //keep finals? apply elsewhere as well..?
		final LaguerreSolver solver = new LaguerreSolver(0.5*ROUND_TOLERANCE);
		Complex[] eigenvalues = solver.solveAllComplex(coefficients, 0);
		
		for(int i=0; i < eigenvalues.length; i++) {	
//			eigenvalues[i] = new Complex(Matrix.roundPolyTerm(eigenvalues[i].getReal()),   //neutralizing
//										 Matrix.roundPolyTerm(eigenvalues[i].getImaginary()));
			
			eigenvalues[i] = Matrix.roundComplexEntry(eigenvalues[i]);      //THIS IS CAUSING EIGENVECTOR ISSUE
		}
		Arrays.sort(eigenvalues, new Comparator<Complex>() {   //sorting the eigenvalues the way I would sort them
			@Override
			public int compare(Complex c1, Complex c2) {
				if(c1.getImaginary() == 0 && c2.getImaginary() != 0) {
					return -1;
				} else if(c1.getImaginary() != 0 && c2.getImaginary() == 0) {       //THIS BRACKET STYLING???????????
					return 1;
				} else {
					if(c1.getReal()-c2.getReal() < -LAGUERRE_TOLERANCE)
						return -1;
					else if(Precision.equals(c1.getReal(), c2.getReal(), LAGUERRE_TOLERANCE))
						if(c1.getImaginary() < c2.getImaginary())
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
			
			//REPLACED ROUND_TOLERANCE WITH 2*POLY_TOLERANCE (ROUND_TOLERANCE FINALLY FAILED.....)
			
			if(Complex.equals(eigenvalues[i], eigenvalues[i-1], 2*Matrix.LAGUERRE_TOLERANCE)) {   //grouping algebraic multiplicities
				
				//System.out.println("TRUE:" + eigenvalues[i].toString() + " equals " + eigenvalues[i-1].toString()); 
				
				//eigenvalues[i-1] = Matrix.roundComplexEntry(eigenvalues[i-1]);
				eigenvalues[i] = eigenvalues[i-1];      //USE CONTAINS INSTEAD?????????????
				distinctEigenvalues.replace(eigenvalues[i], (++multiplicity)-1, multiplicity);   //better (easier to understand) way...?
			} else {
				multiplicity = 1;
				distinctEigenvalues.put(eigenvalues[i],  multiplicity);
			}
		}
		return distinctEigenvalues;
	}
	public LinkedHashSet<Complex[]> eigenvectors(Complex eigenvalue) throws RuntimeException {    //CLEAN UP & TEST
		
		//MAYBE JUST IMPLEMENT THIS METHOD IN ComplexMatrix & CALL FROM HERE (same w/ eigenvalues...?)
		
		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Eigenvectors not defined.");      //IS THIS EVEN TRUE???????
		
		ComplexMatrix complexCopy = new ComplexMatrix(this.entries);	
		ComplexMatrix diffMatrix = complexCopy.minus(new ComplexMatrix(this.numCols).scalarMultiply(eigenvalue));
		diffMatrix = diffMatrix.toRREF();
		Complex[][] diffEntries = diffMatrix.getEntriesCopy();    //USE THIS CORRECTLY!!!!!!!!!!!!!!!!!
		Complex[][] nullSpaceEntries = new Complex[this.numRows][this.numCols];
		
		for(int i=0; i < this.numRows; i++) {         //THINK ABOUT MAKING THIS MORE CLEAR TO READ
			for(int j=0; j < this.numCols; j++) {
				if(!diffEntries[i][j].equals(Complex.ZERO)) {
					nullSpaceEntries[j] = diffEntries[i];
					
					for(int k=0; k < nullSpaceEntries[j].length; k++) {
						if(k == j) {
							nullSpaceEntries[j][k] = Complex.ZERO;
						} else {
							nullSpaceEntries[j][k] = nullSpaceEntries[j][k].multiply(-1.0);
						}
					}
					break;
				}
			}
		}		
		for(int i=0; i < nullSpaceEntries.length; i++) {
			if(nullSpaceEntries[i][0] == null) {                                          //THIS LINE FIXED IT
				Complex[] rowVector = new Complex[this.numCols];
				
				for(int j=0; j < rowVector.length; j++) {
					if(j == i) {
						rowVector[j] = Complex.ONE;
					} else {
						rowVector[j] = Complex.ZERO;
					}
				}			
				nullSpaceEntries[i] = rowVector;
			}
		}
		ComplexMatrix nullSpaceMatrix = new ComplexMatrix(nullSpaceEntries);
		nullSpaceMatrix = nullSpaceMatrix.transpose();
		nullSpaceMatrix.roundEntries();
		Complex[][] eigenSpace = nullSpaceMatrix.getEntriesCopy();
		LinkedHashSet<Complex[]> eigenvectors = new LinkedHashSet<Complex[]>(nullSpaceMatrix.rank());
		
		for(int i=0; i < eigenSpace.length; i++) {
			for(int j=0; j < eigenSpace[i].length; j++) {
				if(!eigenSpace[i][j].equals(Complex.ZERO)) {
					eigenvectors.add(eigenSpace[i]);
					break;
				}
			}
		}
		return eigenvectors;
	}
	
	public Matrix toREF() {   //row echelon form

		double[][] originalEntries = this.getEntriesCopy();
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
	private boolean clearColumnDown(int colIndex, int startRowIndex) {    //helper method (toREF)... MAKE VOID???

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

		double[][] originalEntries = this.getEntriesCopy();
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
	public ComplexMatrix toD() throws RuntimeException {     //diagonal rep consisting of eigenvalues (AEAP...?)

		if(this.numRows != this.numCols) 
			throw new RuntimeException("Matrix isn't square. Diagonal Representation not defined.");

		LinkedHashMap<Complex,Integer> eigenvalues = this.eigenvalues();
		Set<Complex> distinctEigenvalues = eigenvalues.keySet();
		Complex[][] diagonalEntries = new Complex[this.numRows][this.numCols];

		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				diagonalEntries[i][j] = Complex.ZERO;      //Better way to do this...?
		}
		int eigenvalueCount = 0;

		for(Complex eigenvalue: distinctEigenvalues) {
			for(int i=0; i < eigenvalues.get(eigenvalue); i++)
				diagonalEntries[eigenvalueCount][eigenvalueCount++] = eigenvalue;
		}		
		return new ComplexMatrix(diagonalEntries);
	}
	public ComplexMatrix toP() throws RuntimeException {     //CURRENTLY IN "POST-CONSTRUCTION" PHASE
		
		if(this.numRows != this.numCols)
			throw new RuntimeException("Matrix isn't square. Diagonal Representation not defined.");
		
		Complex[][] basisEntries = new Complex[this.numRows][this.numCols];
		int eigenvectorCount = 0;
		Set<Complex> eigenvalues = this.eigenvalues().keySet();
		
		for(Complex eigenvalue: eigenvalues) {
			LinkedHashSet<Complex[]> associatedEigenvectors = this.eigenvectors(eigenvalue);
			for(Complex[] associatedEigenvector: associatedEigenvectors)
				basisEntries[eigenvectorCount++] = associatedEigenvector;           //EROOR HERE (ArrayIndexOutOfBounds)
		}
		if(eigenvectorCount < this.numCols)     //DOUBLE CHECK ON THIS CONDITION
			throw new RuntimeException("Algebraic and geometric multiplicities do not match. "
										+ "Diagonal Representation not defined.");
		
		return (new ComplexMatrix(basisEntries)).transpose();
	}

	/*
	 * TODO: MAKE reduceToREF & reduceToRREF (void) and/or Singular Value Decomposition (SVD)
	 */
	
	/*
	 * //////////////// CONSTRUCTION ZONE BELOW ///////////////////////////////////////////////////////////////////////
	 * 
	 * (none currently)
	 */

	public void setEntries(double[][] entries) {       //START USING THIS IN TESTING & APPLICATION!!!!!!!!!!!!!!!!
		this.entries = entries;
	}
	
	public void setEntry(int i, int j, double entry) {
		this.entries[i][j] = entry;
	}
	
	private static Complex roundComplexEntry(Complex entry) {    //MOVE THIS TO ComplexMatrix (...?) & REPLACE FOR (many) INSTANCES
		double nearestRealWholeNum = Math.round(entry.getReal());
		double nearestImagWholeNum = Math.round(entry.getImaginary());
		Complex roundedEntry = new Complex(nearestRealWholeNum, nearestImagWholeNum);
		
//		if(Complex.equals(entry, roundedEntry, Matrix.POLY_TOLERANCE)) {
//			return roundedEntry;
//		}
		
		if(entry.getReal()-nearestRealWholeNum < Matrix.LAGUERRE_TOLERANCE &&               //10*POLY_TOLERNACE (.0001) WORKS
		   entry.getImaginary()-nearestImagWholeNum < Matrix.LAGUERRE_TOLERANCE) {
			return roundedEntry;
		}
		
		return entry;
	}
	
	
	/*
	 * ///////////// "BEHIND-THE-SCENES METHODS" BELOW ////////////////////////////////////////////////////////////////
	 */

	public static double roundPolyTerm(double entry) {
		double nearestWholeNum = Math.round(entry);
		if(Math.abs(nearestWholeNum-entry) < Matrix.LAGUERRE_TOLERANCE)    //MAKE THIS DISTINCT (NAMING)
			return nearestWholeNum;
		return entry;
	}
	private static double roundEntry(double entry) {
		double nearestWholeNum = Math.round(entry);
		if(Math.abs(nearestWholeNum-entry) < Matrix.ROUND_TOLERANCE)     //MAKE THIS DISTINCT (NAMING)
			return nearestWholeNum;
		return entry;
	}
	
	private void roundEntries() {
		for(int i=0; i < this.numRows; i++) {
			for(int j=0; j < this.numCols; j++)
				this.entries[i][j] = Matrix.roundEntry(this.entries[i][j]);
		}
	}

	//TODO: MAKE SOME SETTERS...?

	public double[][] getEntriesCopy() {
		double[][] entriesCopy = new double[this.entries.length][this.entries[0].length];
		for(int i=0; i < this.entries.length; i++) {
			for(int j=0; j < this.entries[0].length; j++)     //ACTUALLLY saving the original entries
				entriesCopy[i][j] = this.entries[i][j];
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

	public static double evaluatePolynomial(double[] coefficients, double value) {
		double result = 0;
		for(int i=0; i < coefficients.length; i++) {
			result += coefficients[i]*Math.pow(value,i);
		}
		return result;
	}

	public static void printArray(double[] entries) {
		String result = "[";
		for(int i=0; i < entries.length; i++) {
			result += entries[i] + ", ";
		}
		result += "]";
		System.out.println(result);
	}
	public static String arrayString(double[] entries) {
		String result = "[";
		for(int i=0; i < entries.length; i++) {
			result += String.format("% 8.3f", entries[i]) + ", ";
		}
		result += "]";
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
	public static String complexArrayString(Complex[] roots) {
		String result = "[";
		for(int i=0; i < roots.length; i++) {
			result += String.format("% 8.3f", roots[i].getReal()) + " + " 
					+ String.format("% 8.3f", roots[i].getImaginary()) + "i, ";
		}
		result += "]";
		return result;
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
	public static void print2DComplexArray(Complex[][] entries) {
		String result = "[";
		for(int i=0; i < entries.length; i++) {
			result += "[";
			for(int j=0; j < entries[0].length; j++)
				result += entries[i][j].toString() + ", ";
			result += "]";
		}
		result += "]";
		System.out.println(result);
	}

	@Override
	public String toString() {                       //basic toString
		this.roundEntries();
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
		this.roundEntries();
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