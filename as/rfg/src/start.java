import java.math.BigDecimal;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;

public class start {

	/**
	 * @param args
	 */

	static public double Zaokragl_1(double z) {
		z = new BigDecimal(z).setScale(2, BigDecimal.ROUND_HALF_UP)
				.doubleValue();
		return z;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double[][] matrix = { { 1, 0, 0, 1, 0, 0, 0, 0, 0 },
				{ 1, 0, 1, 0, 0, 0, 0, 0, 0 }, { 1, 1, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 1, 1, 0, 1, 0, 0, 0, 0 }, { 0, 1, 1, 2, 0, 0, 0, 0, 0 },
				{ 0, 1, 0, 0, 1, 0, 0, 0, 0 }, { 0, 1, 0, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 1, 1, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0, 0, 0, 1 },
				{ 0, 0, 0, 0, 0, 1, 1, 1, 0 }, { 0, 0, 0, 0, 0, 0, 1, 1, 1 },
				{ 0, 0, 0, 0, 0, 0, 0, 1, 1 } };

		// System.out.print(U);

		double[] suma = new double[9];
		double[] srednia = new double[9];

		System.out.print("Œrednie wektrów : \n");
		for (int i = 0; i < 9; i++) {
			suma[i] = 0;
			srednia[i] = 0;
			for (int j = 0; j < 12; j++) {
				suma[i] += matrix[j][i];
			}
			srednia[i] = suma[i] / 12.0;
			System.out.print("Œrednia[" + i + "] : " + srednia[i] + "\n");
		}
		System.out.print("\n \n");

		System.out.print("Po usuniêciu œredniej : \n");
		double[][] matrixNew = new double[12][9];
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 9; j++) {
				matrixNew[i][j] += matrix[i][j] - srednia[j];
				System.out.print(Zaokragl_1(matrixNew[i][j]) + ", ");
			}
			;
			System.out.print("\n");
		}
		System.out.print("\n \n");

		// double[][] matrixNew = new double[12][9];
		// for(int i=0; i < 12; i++)
		// {
		// for(int j=0; j<9; j++)
		// {
		// matrixNew[i][j] += matrix[i][j] - srednia[j];
		// System.out.print( matrixNew[i][j] + " = " + matrix[i][j] + " - " +
		// srednia[j] + "\n");
		// };
		// }
		// double[] wektor = new double[12];

		// for(int i=0; i < 12; i++)
		// {
		// wektor[i] = 0;
		// for(int j=0; j<9; j++)
		// {
		// wektor[i] = matrixNew[i][j];
		// }
		// }

		// double[] wektor = new double[12];
		// wektor[0] = 0;
		// for(int i =0; i<12; i++)
		// {
		// wektor[0] += matrixNew[i][0] * matrixNew[i][1];
		// }
		// System.out.print("C1/C2 = " + wektor[0]);

		double dll = 0;
		for (int i = 0; i < 12; i++) {
			dll += Math.pow(matrixNew[i][0], 2);
		}
		double a = Math.sqrt(dll);
		System.out.print("d1: " + Math.sqrt(dll));

		dll = 0;
		for (int i = 0; i < 12; i++) {
			dll += Math.pow(matrixNew[i][1], 2);
		}
		double b = Math.sqrt(dll);
		System.out.print("d2: " + Math.sqrt(dll) + "\n");
		// System.out.print("il: " + a * b + " = " + wektor[0] / (a*b) );

		double[][] wektor = new double[9][9];
		for (int j = 0; j < 8; j++) {
			for (int k = j + 1; k < 9; k++) {
				for (int i = 0; i < 12; i++) {
					wektor[j][k] += matrixNew[i][j] * matrixNew[i][k];
				}
			}
		}
		WypiszMatrixa(wektor, "Wektory", 9, 9);

		double[] dl = new double[9];
		for (int j = 0; j < 9; j++) {
			for (int i = 0; i < 12; i++) {
				dl[j] += Math.pow(matrixNew[i][j], 2);
			}
			dl[j] = Math.sqrt(dl[j]);
		}

		for (int j = 0; j < 8; j++) {
			for (int k = j + 1; k < 9; k++) {
				wektor[j][k] = wektor[j][k] / (dl[k] * dl[j]);
			}
		}

		double[][] przedSVD = wektor;
		WypiszMatrixa(wektor, "RELACJA MIÊDZY DOKUMENTAMI PRZED SVD", 9, 9);

		// System.out.print("d2: " + Math.sqrt(dl) + "\n");
		// System.out.print("il: " + a * b + " = " + wektor[0] / (a * b));

		/*
		 * 
		 * ZADANIE 2 PO SVD
		 */

		DenseDoubleMatrix2D ok = new DenseDoubleMatrix2D(matrix);
		SingularValueDecomposition c = new SingularValueDecomposition(ok);
		DoubleMatrix2D U = c.getU();
		DoubleMatrix2D S = c.getS();
		DoubleMatrix2D V = c.getV();

		double[][] matrixU = toMatrix(U);
		double[][] matrixV = toMatrix(V);
		double[][] matrixS = toMatrix(S);

		double[][] newMatrixU = new double[12][2];
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 2; j++) {
				newMatrixU[i][j] = matrixU[i][j];
			}
		}
		
		double[][] newMatrixV = new double[2][9];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 9; j++) {
				newMatrixV[i][j] = matrixV[i][j];
			}
		}
		
		double[][] newMatrixS = new double[2][2];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				newMatrixS[i][j] = matrixS[i][j];
			}
		}
	}

	public static double[][] toMatrix(DoubleMatrix2D V) {
		double[][] mat = new double[V.rows()][V.columns()];
		for (int i = 0; i < V.rows(); i++) {
			for (int j = 0; j < V.columns(); j++) {
				mat[i][j] = V.get(i, j);
			}
		}
		return mat;
	}

	public static int[][] transponuj(int x[][]) {
		int[][] y;

		y = new int[x.length][x.length];
		for (int row = 0; row < x.length; row++) {
			for (int column = 0; column < x[row].length; column++) {
				y[column][row] = x[row][column];
			}
		}
		return y;
	}

	private static void BezWypiszMatrixa(double[][] wektor, String tytul,
			int wiersze, int kolumny) {
		System.out.print(tytul + " : \n");
		for (int i = 0; i < wiersze; i++) {
			for (int j = 0; j < kolumny; j++) {
				System.out.print(wektor[i][j] + ", ");
			}
			;
			System.out.print("\n");
		}
		System.out.print("\n \n");

	}

	private static void WypiszMatrixa(double[][] wektor, String tytul,
			int wiersze, int kolumny) {
		System.out.print(tytul + " : \n");
		for (int i = 0; i < wiersze; i++) {
			for (int j = 0; j < kolumny; j++) {
				System.out.print(Zaokragl_1(wektor[i][j]) + ", ");
			}
			;
			System.out.print("\n");
		}
		System.out.print("\n \n");

	}
}
