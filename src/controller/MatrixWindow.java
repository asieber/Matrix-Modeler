package controller;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import org.apache.commons.math3.complex.Complex;

import model.Matrix;

/*
 * 
 * A class that displays the following information about the matrix, A, entered by the user:
 * 
 * A
 * tr(A)
 * det(A)
 * rank(A)
 * AT
 * AREF
 * ARREF
 * AInv
 * range(A)
 * nullspace(A)
 * charPoly(A)
 * eigenvalues(A) w/ associated multiplicities
 * eigenvectors(A) w/ associated eigenvalues
 * PDPINV
 * 
 */

public class MatrixWindow extends JFrame 
{
	private double[][] AEntries 			= null;	               //THINK OF MOVING THESE VARIABLES ELSEWHERE
	private Matrix A 						= null;
	private final static int UNINITIALIZED  = -1;
	private int numRows 					= UNINITIALIZED;
	private int numCols 					= UNINITIALIZED;
	private final static String INDENT 		= "                    ";
	
	public MatrixWindow() 
	{
		layoutFrame();
	}
	
	private void layoutFrame()      
	{
		this.setLocation(0,0);
		this.setSize(1360, 730);
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				
		JPanel view 				= new JPanel();
		JPanel topPanel 			= new JPanel();
		JPanel topLeftPanel 		= new JPanel();
		JPanel topRightPanel 		= new JPanel();
		JPanel topRightLeftPanel 	= new JPanel();
		JPanel topRightRightPanel 	= new JPanel();
		JPanel separatorPanel 		= new JPanel();
		
		class MatrixPanel extends JPanel {    //Essential for organization/style
			
			JTextArea rightText;
			
			public MatrixPanel(String name) {
				
				JPanel leftPanel 	= new JPanel();
				JPanel rightPanel 	= new JPanel();
				rightText 			= new JTextArea("A " + name);
				
				this.setLayout(new GridLayout(1,2));
				leftPanel.setLayout(new GridLayout(1,1));
				rightPanel.setLayout(new GridLayout(1,1));
				
				this.setPreferredSize(rightText.getSize());
				
				JLabel leftText = new JLabel(MatrixWindow.INDENT + "A " + name + " = ........................"
																			+ "...................");
				
				leftText.setFont(new Font("Courier", Font.BOLD, 12));
				leftText.setHorizontalAlignment(SwingConstants.LEFT);
				leftText.setVerticalAlignment(SwingConstants.TOP);
				leftPanel.add(leftText);
				
				rightText.setFont(new Font("Courier", 0, 10));
				rightPanel.add(rightText);
				this.add(leftPanel);
				this.add(rightPanel);
			}
			
			public void setRightText(String text) {
				rightText.setText(text);
			}
		}
		
		MatrixPanel APanel 				= new MatrixPanel("");
		MatrixPanel ATracePanel 		= new MatrixPanel("trace");
		MatrixPanel ADetPanel 			= new MatrixPanel("determinant");
		MatrixPanel ARankPanel 			= new MatrixPanel("rank");
		MatrixPanel ATransPanel 		= new MatrixPanel("transpose");
		MatrixPanel AREFPanel 			= new MatrixPanel("row echelon form");
		MatrixPanel ARREFPanel 			= new MatrixPanel("reduced row echelon form");
		MatrixPanel AInvPanel 			= new MatrixPanel("inverse");
		MatrixPanel ARangePanel 		= new MatrixPanel("range");
		MatrixPanel ANullPanel 			= new MatrixPanel("nullspace");
		MatrixPanel ACharPolyPanel 		= new MatrixPanel("characteristic polynomial");
		MatrixPanel AEigenvaluesPanel 	= new MatrixPanel("eigenvalues");
		MatrixPanel AEigenvectorsPanel 	= new MatrixPanel("eigenvectors");
		MatrixPanel ADiagonalPanel 		= new MatrixPanel("diagonalized");
		
		view.setLayout(new GridLayout(16,1));
		topPanel.setLayout(new GridLayout(1,2));
		topLeftPanel.setLayout(new GridLayout(2,2));
		topRightPanel.setLayout(new GridLayout(1,2));
		topRightLeftPanel.setLayout(new GridLayout(1,1));
		topRightRightPanel.setLayout(new GridLayout(1,1));
		separatorPanel.setLayout(new GridLayout(1,1));

		JLabel entriesTitle = new JLabel(MatrixWindow.INDENT + "Enter the entries (Must press enter after each one): ");
		entriesTitle.setFont(new Font("Courier", Font.BOLD, 12));
		topRightLeftPanel.add(entriesTitle);
		JLabel rowTitle = new JLabel(MatrixWindow.INDENT + "Enter the number of rows (Must press enter after): ");
		rowTitle.setFont(new Font("Courier", Font.BOLD, 12));
		topLeftPanel.add(rowTitle);
		
		class EntryTextListener implements ActionListener
		{
			JTextField[][] entryFields;
			int i, j;
			
			EntryTextListener(JTextField[][] entryFields, int i, int j) 
			{
				this.entryFields = entryFields;
				this.i = i;
				this.j = j;
			}
			@Override
			public void actionPerformed(ActionEvent e) 
			{
				String entryString = this.entryFields[this.i][this.j].getText();
				
				if(!entryString.equals("") && entryString != null) {
					AEntries[this.i][this.j] = Double.parseDouble(entryString);
				}
				
				A = new Matrix(AEntries);
				APanel.setRightText(A.toString());
				
				try {
					double ATrace = (new Matrix(A.getEntriesCopy())).trace();
					ATracePanel.setRightText("" + ATrace);
				} catch(RuntimeException rte) {
					ATracePanel.setRightText(rte.getMessage());
				}
				
				try {
					double ADet = (new Matrix(A.getEntriesCopy())).det();
					ADetPanel.setRightText("" + ADet);
				} catch(RuntimeException rte) {
					ADetPanel.setRightText(rte.getMessage());
				}
				
				int ARank = (new Matrix(A.getEntriesCopy())).rank();
				ARankPanel.setRightText("" + ARank);
				
				Matrix ATrans = (new Matrix(A.getEntriesCopy())).transpose();
				ATransPanel.setRightText(ATrans.toString());
				
				Matrix AREF = (new Matrix(A.getEntriesCopy())).toREF();
				AREFPanel.setRightText(AREF.toString());
				
				Matrix ARREF = (new Matrix(A.getEntriesCopy())).toRREF();     //NEED Copy
				ARREFPanel.setRightText(ARREF.toString());
				
				try {
					Matrix AInv = (new Matrix(A.getEntriesCopy())).inverse();
					AInvPanel.setRightText(AInv.toString());
				} catch(RuntimeException rte) {
					AInvPanel.setRightText(rte.getMessage());
				}
				
				LinkedHashSet<double[]> range = (new Matrix(A.getEntriesCopy())).range();
				String rangeBasis = "";
				for(double[] basisVector: range) {
					rangeBasis += Matrix.arrayString(basisVector) + "\n";
				}
				ARangePanel.setRightText(rangeBasis);
			
				//CONSIDER IMPLEMENTING Matrix.nullspace() (SINCE JUST DOUBLES, NOT COMPLEX)
				
				LinkedHashSet<double[]> ANullspace = (new Matrix(A.getEntriesCopy())).nullspace();
				String nullspaceString = "";
				for(double[] vector: ANullspace) {
					nullspaceString += Matrix.arrayString(vector) + "\n";
				}
				if(nullspaceString.equals(""))
					ANullPanel.setRightText("Nullspace is empty.");
				else
					ANullPanel.setRightText(nullspaceString);
				
				try {
					double[] ACharPoly = (new Matrix(A.getEntriesCopy())).characteristicPolynomial();
					String charPoly = "";
					for(int i=ACharPoly.length-1; i > -1; i--) {
						charPoly += "" + ACharPoly[i] + "x^" + i + " + ";
					}
					ACharPolyPanel.setRightText(charPoly);
				} catch(RuntimeException rte) {
					ACharPolyPanel.setRightText(rte.getMessage());
				}
				
				try {
					LinkedHashMap<Complex,Integer> AEigenvalues = (new Matrix(A.getEntriesCopy())).eigenvalues();
					Set<Complex> eigenvalues = AEigenvalues.keySet();
					String eigenvaluesString = "";
					for(Complex eigenvalue: eigenvalues) {
						eigenvaluesString += String.format("% 8.3f", eigenvalue.getReal()) + " + " 
										   + String.format("% 8.3f", eigenvalue.getImaginary()) + "i"   //String.format
										   + " (multiplicity: " + AEigenvalues.get(eigenvalue) + ")\n";
					}
					AEigenvaluesPanel.setRightText(eigenvaluesString);
				} catch(RuntimeException rte) {
					AEigenvaluesPanel.setRightText(rte.getMessage());
				}
								
				try {
										
					LinkedHashMap<Complex,Integer> AEigenvalues = (new Matrix(A.getEntriesCopy())).eigenvalues();
					Set<Complex> eigenvalues = AEigenvalues.keySet();
					String eigenvectorsString = "";
					for(Complex eigenvalue: eigenvalues) {
						LinkedHashSet<Complex[]> AEigenvectors = (new Matrix(A.getEntriesCopy())).eigenvectors(eigenvalue);
						for(Complex[] eigenvector: AEigenvectors) {
							eigenvectorsString += Matrix.complexArrayString(eigenvector) + " (eigenvalue: "
												+ String.format("% 8.3f", eigenvalue.getReal()) + " + " 
												+ String.format("% 8.3f", eigenvalue.getImaginary()) + "i)\n";
						}
					}
					AEigenvectorsPanel.setRightText(eigenvectorsString);
				} catch(RuntimeException rte) {
					AEigenvectorsPanel.setRightText(rte.getMessage());
				}
				
				try {
					String diagonalizedRep = "";
					diagonalizedRep += (new Matrix(A.getEntriesCopy())).toP().toString() + "\n";
					diagonalizedRep += (new Matrix(A.getEntriesCopy())).toD().toString() + "\n";
					diagonalizedRep += (new Matrix(A.getEntriesCopy())).toP().inverse().toString();
					ADiagonalPanel.setRightText(diagonalizedRep);
				} catch(RuntimeException rte) {
					ADiagonalPanel.setRightText(rte.getMessage());
				}
			}
		}
		
		JTextField numRowsField = new JTextField();
		
		class RowTextListener implements ActionListener 
		{
			@Override
			public void actionPerformed(ActionEvent e) 
			{
				String numRowsString = numRowsField.getText();
				
				if(!numRowsString.equals("") && numRowsString != null) {
					numRows = (int) Double.parseDouble(numRowsString);
				}
				
				if(numRows != MatrixWindow.UNINITIALIZED && numCols != MatrixWindow.UNINITIALIZED) {
					
					AEntries = new double[numRows][numCols];
					A = new Matrix(AEntries);
					
					//initialize new entryFields
					
					JTextField[][] entryFields = new JTextField[numRows][numCols];
					
					for(int i=0; i < numRows; i++) {
						for(int j=0; j < numCols; j++) {
							entryFields[i][j] = new JTextField();    //MOVE IN SOMEWHERE...?
						}
					}
					
					//remove the old TextFields
					
					topRightRightPanel.setLayout(new GridLayout(numRows, numCols));
					topRightRightPanel.removeAll();
					
					//add the new TextFields
					
					for(int i=0; i < numRows; i++) {
						for(int j=0; j < numCols; j++) {
							topRightRightPanel.add(entryFields[i][j]);
						}
					}
					
					//add listeners
					
					for(int i=0; i < numRows; i++) {
						for(int j=0; j < numCols; j++) {
							entryFields[i][j].addActionListener(new EntryTextListener(entryFields,i,j));
						}
					}
				}
				
				if(A != null) {
					APanel.setRightText(A.toString());
				}
			}
		}
				
		numRowsField.addActionListener(new RowTextListener());
		topLeftPanel.add(numRowsField);
		JLabel colTitle = new JLabel(MatrixWindow.INDENT + "Enter the number of cols (Must press enter after): ");   //CHNAGE BACK TO "columns"...
		colTitle.setFont(new Font("Courier", Font.BOLD, 12));
		topLeftPanel.add(colTitle);
		JTextField numColsField = new JTextField();    //add Listener...?
		
		class ColTextListener implements ActionListener 
		{
			@Override
			public void actionPerformed(ActionEvent e) 
			{
				String numColsString = numColsField.getText();
				
				if(!numColsString.equals("") && numColsString != null) {
					numCols = (int) Double.parseDouble(numColsString);
				}
				
				if(numRows != MatrixWindow.UNINITIALIZED && numCols != MatrixWindow.UNINITIALIZED) {
					
					AEntries = new double[numRows][numCols];
					A = new Matrix(AEntries);
					
					//initialize new entryFields
					
					JTextField[][] entryFields = new JTextField[numRows][numCols];
					
					for(int i=0; i < numRows; i++) {
						for(int j=0; j < numCols; j++) {
							entryFields[i][j] = new JTextField();    //MOVE IN SOMEWHERE...?
						}
					}
					
					//remove the old TextFields
					
					topRightRightPanel.setLayout(new GridLayout(numRows, numCols));					
					topRightRightPanel.removeAll();
					
					//add tihe new TextFields
					
					for(int i=0; i < numRows; i++) {
						for(int j=0; j < numCols; j++) {
							topRightRightPanel.add(entryFields[i][j]);
						}
					}
					
					//add Listeners
					
					for(int i=0; i < numRows; i++) {
						for(int j=0; j < numCols; j++) {
							entryFields[i][j].addActionListener(new EntryTextListener(entryFields,i,j));
						}
					}
				}
				
				if(A != null) {
					APanel.setRightText(A.toString());
				}
			}
		}
				
		numColsField.addActionListener(new ColTextListener());
		topLeftPanel.add(numColsField);
		
		topPanel.add(topLeftPanel);
		topRightPanel.add(topRightLeftPanel);
		topRightPanel.add(topRightRightPanel);
		topPanel.add(topRightPanel);
		
		view.add(topPanel);
		view.add(separatorPanel);		
		view.add(APanel);
		view.add(ATracePanel);
		view.add(ADetPanel);
		view.add(ARankPanel);
		view.add(ATransPanel);
		view.add(AREFPanel);
		view.add(ARREFPanel);
		view.add(AInvPanel);
		view.add(ARangePanel);
		view.add(ANullPanel);
		view.add(ACharPolyPanel);
		view.add(AEigenvaluesPanel);
		view.add(AEigenvectorsPanel);
		view.add(ADiagonalPanel);
		
		JScrollPane scrollpane = new JScrollPane(view, 
												 JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, 
												 JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		this.add(scrollpane);
		this.setVisible(true);
	}
	
	public static void main(String[] args) 
	{
		new MatrixWindow();
	}

}