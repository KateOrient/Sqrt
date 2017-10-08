package Sqrt;

import java.io.*;
import java.util.StringTokenizer;

public class Sqrt{
    private double[][] matrix;
    private double[] f;
    private int n;
    private double determinant;

    public Sqrt (String fileName) throws IOException{
        loadFromFile(fileName);
    }

    public void loadFromFile (String fileName) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        n = Integer.parseInt(reader.readLine());
        matrix = new double[n][n];
        f = new double[n];
        String s;
        for (int i = 0; i < n; i++){
            s = reader.readLine();
            StringTokenizer st = new StringTokenizer(s);
            for (int j = 0; j < n; j++){
                matrix[i][j] = Double.parseDouble(st.nextToken());
            }
            f[i] = Double.parseDouble(st.nextToken());
        }
        reader.close();
    }

    public void printMatrix (){
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                System.out.format("%10.05f ", matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    public void printF (){
        printVector(f, n);
    }

    //транспонирование матрицы
    public double[][] transpose (){
        double[][] tr = new double[n][n];
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                tr[j][i] = matrix[i][j];
            }
        }
        return tr;
    }

    //домножение матрицы А и неоднородности на Ат(транспонированная)
    public void makeSymetric (){
        double[][] sym = new double[n][n];
        double[][] tr = transpose();
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                for (int k = 0; k < n; k++){
                    sym[i][j] += tr[i][k] * matrix[k][j];
                }
            }
        }
        double[] atf = new double[n];
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                atf[i] += tr[i][j] * f[j];
            }
        }
        matrix = sym;
        f = atf;
    }

    public void printVector (double[] v, int r){
        System.out.print("(");
        for (int i = 0; i < r; i++){
            System.out.format("%10.05f", v[i]);
            if (i != r - 1){
                System.out.print(", ");
            }
        }
        System.out.println(")");
    }

    //решение системы методом квадратного корня
    public double[] solve (){
        double[] x = new double[n];
        double[] y = new double[n];
        double[][] s = new double[n][n];
        double[] d = new double[n];
        d[0] = matrix[0][0] / Math.abs(matrix[0][0]);

        //нахождения треугольной матрицы S
        s[0][0] = Math.sqrt(matrix[0][0] * d[0]);
        for (int j = 1; j < n; j++){
            s[0][j] = matrix[0][j] / s[0][0] * d[0];
        }
        for (int i = 1; i < n; i++){
            s[i][i] = matrix[i][i];
            for (int k = 0; k < i; k++){
                s[i][i] -= Math.pow(s[k][i], 2);
            }
            d[i]=s[i][i] / Math.abs(s[i][i]);
            s[i][i] = Math.sqrt(s[i][i]*d[i]);

            if (i != n - 1){
                for (int j = i + 1; j < n; j++){
                    s[i][j] = matrix[i][j];
                    for (int k = 0; k < i; k++){
                        s[i][j] -= s[k][i] * s[k][j]*d[k];
                    }
                    s[i][j] /= s[i][i]*d[i];
                }
            }
        }
        System.out.println("Треугольная матрица S:");
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                System.out.format("%10.05f ", s[i][j]);
            }
            System.out.println();
        }
        System.out.println();
        //нахождение вектора y
        y[0] = f[0] / s[0][0];
        for (int i = 1; i < n; i++){
            y[i] = f[i];
            for (int k = 0; k < i; k++){
                y[i] -= s[k][i] * y[k];
            }
            y[i] /= s[i][i];
        }

        System.out.print("Вектор y:");
        printVector(y, n);

        //нахождение вектора х
        x[n - 1] = y[n - 1] / (s[n - 1][n - 1]);
        for (int i = n - 2; i >= 0; i--){
            x[i] = y[i];
            for (int k = i + 1; k < n; k++){
                x[i] -= s[i][k] * x[k]*d[i];
            }
            x[i] /= s[i][i]*d[i];
        }

        System.out.print("Вектор x:");
        printVector(x, n);

        //высчитывание определтеля матрицы
        determinant = 1;
        for(int i =0;i<n;i++){
            determinant*=s[i][i];
        }

        return x;
    }


    //проверка вектора рещения
    public double[] checkSol (double[] sol){
        double[] discrepancy = new double[n];
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                discrepancy[i] += matrix[i][j] * sol[j];
            }
            discrepancy[i] -= f[i];
        }
        return discrepancy;
    }


    public double getDeterminant (){
        return determinant;
    }

    public int getN (){
        return n;
    }
}
