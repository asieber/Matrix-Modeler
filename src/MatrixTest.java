import static org.junit.Assert.*;

import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;
import org.junit.Test;

public class MatrixTest {

	/*
	@Test
	public void testMatrixToStringLong() {
		Matrix A = new Matrix(new double[][] {{1,2,3},
											  {4,5,6},
											  {7,8,9}});
		System.out.println("A toStringLong: ");
		System.out.println(A.toStringLong());
	}
	 */

	/*
	@Test
	public void testMatrixSquared() {
		Matrix A = new Matrix(new double[][] {{0,2,3},
											  {0,5,6},
											  {0,8,9}});
		Matrix AA = A.rightMultiply(A);
		System.out.println("AA toString: ");
		System.out.println(AA.toString());
	}
	 */

	/*
	@Test
	public void testMatrixREF() {
		Matrix A = new Matrix(new double[][] {{3,4,-3},
			  								  {2,5,5},
			  								  {1,-3,1}});
		Matrix refA = A.toREF();
		System.out.println("A REF toString: ");
		System.out.println(refA.toString());
		System.out.println("rank = " + refA.rank());
	}

	@Test
	public void testMatrixRREF() {
		Matrix A = new Matrix(new double[][] {{3,4,-3},
			  								  {2,5,5},
			  								  {1,-3,1}});    //-2,-3,1
		Matrix rrefA = A.toRREF();
		System.out.println("A RREF toString: ");
		System.out.println(rrefA.toString());
		System.out.println("rank = " + rrefA.rank());
	}
	 */

	/*
	@Test
	public void testMatrixTranspose() {
		Matrix A = new Matrix(new double[][] {{-7,-6,-12,-33},
			  {5,5,7,24},
			  {1,0,4,5}});
		Matrix AT = A.transpose();
		System.out.println("AT toString: ");
		System.out.println(AT.toString());
		System.out.println("A toString: ");
		System.out.println(A.toString());
	}
	 */

	/*
	@Test
	public void testMatrixMult() {
		Matrix A = new Matrix(new double[][] {{-7,-6,-12,-33},
			  {5,5,7,24},
			  {1,0,4,5}});
		Matrix AT = A.transpose();
		Matrix AAT = A.rightMultiply(AT);
		System.out.println("AA toString: ");
		System.out.println(AAT.toString());
		System.out.println("A toString: ");
		System.out.println(A.toString());
	}
	 */

	/*
	@Test
	public void testMatrixToString() {
		Matrix A = new Matrix(new double[][] {{-5,10,0,8},
											  {3,-7,1,-4},
											  {-1,0,-2,3},
											  {11,4,-6,0}});
		System.out.println("A toString: ");
		System.out.println(A.toString());
		System.out.println("det(A) = " + A.det() + "\n"); 

		Matrix AINV = A.inverse();
		System.out.println("AINV toString: ");
		System.out.println(AINV.toString());
		System.out.println("det(AINV) = " + AINV.det() + "\n"); 

		Matrix AAINV = A.rightMultiply(AINV);
		System.out.println("AAINV toString: ");
		System.out.println(AAINV.toString());
		System.out.println("det(A)det(AINV) = " + A.det()*AINV.det() + "\n");



		//WRITE OUT WHICH ROW WAS A LIN COMBO OF OTHER & HOW (IF NOT INV)


	}
	 */

	@Test
	public void testQuinticFunction2() {

//		Matrix A = new Matrix(new double[][] {{2,0,0,0,0},    //charPoly: 1,1,-121
//											  {0,1,0,0,0},
//											  {0,0,2,0,0},
//											  {0,0,0,1,0},
//											  {0,0,0,0,1}});

		Matrix A = new Matrix();
		System.out.println("A toString: ");
		System.out.println(A.toString());
		System.out.println("det(A) = " + A.det() + "\n"); 
		System.out.print("charPoly(A) = "); 
		final double[] coefficients = A.characteristicPolynomial();
		Matrix.printArray(coefficients);
		System.out.println("");

		final LaguerreSolver solver = new LaguerreSolver();
		final Complex[] roots = solver.solveAllComplex(coefficients, 0);

		/*
		for (Complex expected : new Complex[] { new Complex(0, -2),
				new Complex(0, 2),
				new Complex(0.5, 0.5 * FastMath.sqrt(3)),
				new Complex(-1, 0),
				new Complex(0.5, -0.5 * FastMath.sqrt(3.0)) }) {
			final double tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
					FastMath.abs(expected.abs() * solver.getRelativeAccuracy()));
			//TestUtils.assertContains(result, expected, tolerance);
		}
		*/
		
		ComplexMatrix D = A.toD();
		System.out.println("D toString: ");
		System.out.println(D.toString());
		
		/*
		for(int i=0; i < roots.length; i++) {
			System.out.println("CharPoly evaluated at x = " + roots[i].getReal() + "\t\tequals\t\t" 
					+ Matrix.evaluate(coefficients, roots[i].getReal()));
		}
		*/
		
	}

	/*
	@Test
	public void testInverseDetAndCharPoly() {

		Matrix A = new Matrix(new double[][] {{-9,7},    //charPoly: 1,1,-121
			{7,8}});

		//Matrix A = new Matrix();
		System.out.println("A toString: ");
		System.out.println(A.toString());
		System.out.println("det(A) = " + A.det() + "\n"); 
		System.out.print("charPoly(A) = "); 
		double[] coefficients = A.characteristicPolynomial();
		Matrix.printArray(coefficients);
		System.out.println("");

		
		//double[] tempCoefficients = coefficients;
		//for(int i=0; i < coefficients.length; i++) {
		//	coefficients[i] = tempCoefficients[coefficients.length-i-1];
		//}
		 

		LaguerreSolver ls = new LaguerreSolver();
		Complex[] roots = ls.solveAllComplex(coefficients, 0);    //should be -11.5 & 10.5
		System.out.print("roots: ");
		Matrix.printComplexArray(roots);
		System.out.println("");

		for(int i=0; i < roots.length; i++) {
			System.out.println("CharPoly evaluated at x = " + roots[i].getReal() + "\t\tequals\t\t" 
					+ Matrix.evaluate(coefficients, roots[i].getReal()));
		}
		*/

		/*
		Matrix AINV = A.inverse();
		System.out.println("AINV toString: ");
		System.out.println(AINV.toString());
		System.out.println("det(AINV) = " + AINV.det() + "\n"); 
		System.out.println("charPoly(AINV) = "); 
		Matrix.printArray(AINV.characteristicPolynomial());
		System.out.println("");

		Matrix AAINV = A.rightMultiply(AINV);
		System.out.println("AAINV toString: ");
		System.out.println(AAINV.toString());
		System.out.println("det(A)det(AINV) = " + A.det()*AINV.det() + "\n");
		System.out.println("charPoly(AAINV) = "); 
		Matrix.printArray(AAINV.characteristicPolynomial());
		System.out.println("");
		

	}
	*/

	/*
	@Test
	public void testMatrixInverse() {
		Matrix A = new Matrix(new double[][] {{3,4,-3},
			  								  {2,5,5},
			  								  {1,-3,1}});
		Matrix AINV = A.inverse();
		System.out.println("AINV toString: ");
		System.out.println(AINV.toString());
	}

	@Test
	public void testMultInverse() {
		Matrix A = new Matrix(new double[][] {{3,4,-3},
			  								  {2,5,5},
			  								  {1,-3,1}});
		Matrix AINV = A.inverse();
		Matrix AINVA = AINV.rightMultiply(A);
		System.out.println("AINVA toString: ");
		System.out.println(AINVA.toString());
	}
	 */

}