import java.util.Arrays;

public class EllipsoidMethod {

    static final double EPSILON = 1e-5;
    static final int MAX_ITER = 1000;

    // Ax <= b
    public static boolean solve(double[][] A, double[] b, double[] x0, double R) {
        int n = x0.length;
        double[][] Q = new double[n][n];

        // Initialize Q = R^2 * I
        for (int i = 0; i < n; i++) Q[i][i] = R * R;

        double[] x = Arrays.copyOf(x0, n);

        for (int iter = 0; iter < MAX_ITER; iter++) {
            boolean allSatisfied = true;

            for (int i = 0; i < A.length; i++) {
                double dot = dot(A[i], x);
                if (dot > b[i] + EPSILON) {
                    allSatisfied = false;

                    double[] a = A[i];
                    double normSq = dot(a, mult(Q, a));

                    double alpha = 1.0 / (n + 1);
                    double[] Qa = mult(Q, a);

                    // Update center
                    for (int j = 0; j < n; j++)
                        x[j] = x[j] - (alpha / Math.sqrt(normSq)) * Qa[j];

                    // Update Q
                    for (int r = 0; r < n; r++) {
                        for (int c = 0; c < n; c++) {
                            Q[r][c] = (n * n / (n * n - 1.0)) * (Q[r][c]
                                    - (2.0 / (n + 1)) * Qa[r] * Qa[c] / normSq);
                        }
                    }
                    break;
                }
            }

            if (allSatisfied) {
                System.out.println("Solution: " + Arrays.toString(x));
                return true;
            }
        }

        System.out.println("No feasible solution found.");
        return false;
    }

    static double dot(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) sum += a[i] * b[i];
        return sum;
    }

    static double[] mult(double[][] M, double[] v) {
        double[] result = new double[v.length];
        for (int i = 0; i < v.length; i++)
            for (int j = 0; j < v.length; j++)
                result[i] += M[i][j] * v[j];
        return result;
    }

    public static void main(String[] args) {
        double[][] A = { {1, 1}, {-1, 2}, {2, 1} };
        double[] b = {4, 2, 5};
        double[] x0 = {0, 0};
        double R = 10;

        solve(A, b, x0, R);
    }
}
