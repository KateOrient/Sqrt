import Sqrt.Sqrt;

import java.io.IOException;

public class Main{
    public static void main (String[] args) throws IOException{
        Sqrt x = new Sqrt("1.txt");
        System.out.println("Исходная матрица системы A:");
        x.printMatrix();
        System.out.println("Вектор неоднородности f:");
        x.printF();
        x.makeSymetric();
        System.out.println("Aт * A:");
        x.printMatrix();
        System.out.println("Ат * f:");
        x.printF();
        double[] sol= x.solve();
        double[] discrepancy = x.checkSol(sol);
        System.out.print("Проверка решения:(");
       for (int i = 0; i < x.getN(); i++){
            System.out.print( discrepancy[i]);
            if (i != x.getN() - 1){
                System.out.print(", ");
            }
        }
        System.out.println(")");
        System.out.format("Определитель матрицы системы: %10.05f", x.getDeterminant());
    }
}
