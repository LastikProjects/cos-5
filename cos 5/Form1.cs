using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace cos_5
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void massiv1()
        {
            double h = Convert.ToDouble(textBox1.Text);
            double hsmall = Convert.ToDouble(textBox2.Text);
            int n = Convert.ToInt32(2 / h) + 1;
            int nbig = Convert.ToInt32(2 / hsmall) + 1;
            double[] x = new double[nbig];
            double[] y = new double[nbig];
            double[] a = new double[n];
            double[] b = new double[n];
            double[] c = new double[n];
            double[] d = new double[n];
            double[] u = new double[n];
            double[] v = new double[n];
            double[] m = new double[n];
            double[] f = new double[n];
            double[] xi = new double[n];
            double[] s = new double[nbig+200];
            xi[0] = -1;
            for(int i=1;i<n;i++)
            {
                xi[i] = xi[i - 1] + h;
            }

            for(int i=0;i<n;i++)
            {
                f[i]= Math.Exp(-xi[i] * xi[i]);//узлы
            }

            x[0] = -1;
            for(int i=1;i<nbig;i++)
            {
                x[i] = x[i - 1] + hsmall;
            }
            
            for(int i=0;i<nbig;i++)//сама функция
            {
                y[i] = Math.Exp(-x[i] * x[i]);
            }
           
            a[0] = 1;
            a[n - 1] = 1;
            b[0] = 0;
            c[0] = 0;
            c[n - 1] = 0;
            d[0] = -2* Math.Exp(-xi[0] * xi[0])+4*x[0]*x[0]* Math.Exp(-xi[0] * xi[0]);
            d[n - 1] = -2 * Math.Exp(-xi[n-1] * xi[n-1]) + 4 * x[n-1] * x[n-1] * Math.Exp(-xi[n-1] * xi[n-1]);
            for (int i = 1; i < n - 1; i++)//матрица для m 
            {
                a[i] = 4;
                b[i] = 1;
                c[i] = 1;
                d[i] = 6 * (f[i + 1] - 2 * f[i] + f[i - 1]) / (h * h);
            }

            u[0] = b[0] / a[0] * (-1);
            v[0] = d[0] / a[0];
            for (int i = 1; i < n; i++)//вспомогательные переменные для m 
            {
                u[i] = (-1) * b[i] / (c[i] * u[i - 1] + a[i]);
                v[i] = (d[i] - c[i] * v[i - 1]) / (c[i] * u[i - 1] + a[i]);
            }

            m[n - 1] = v[n - 1];
            for (int i = n - 2; i > -1; i--)//вычисление m 
            {
                m[i] = u[i] * m[i + 1] + v[i];
            }
            /*//
            double[] masm = new double[nbig];
            double[] masf = new double[nbig];
            for (int i = 1; i < n-1; i++)
            {
                masm[i]=m[i - 1] + 4 * m[i] + m[i + 1];
                masf[i] = 6 * (f[i - 1] + f[i + 1] - 2 * f[i]) / (h * h);
                //masm[i] = masm[i] - masf[i];
            }
            //*///проверка вычислении матрицы
            int k = 0;
            double t;
            for(int i=0;i<n-1;i++)
            {
                for(int j=0;j<h/hsmall;j++)
                {
                    t = hsmall * j / h;
                    s[k] = f[i] * psi1(t) + f[i + 1] * t + m[i] * h * h * psi3(t) + m[i + 1] * h * h * psi4(t);
                    k++;
                }
            }
            t = 1;
            s[nbig - 1] = f[n - 2] * psi1(t) + f[n - 1] * t + m[n - 2] * h * h * psi3(t) + m[n - 1] * h * h * psi4(t);

            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            for (int i = 0; i < nbig; i++)//вывод первой функции 
            {
                chart1.Series[0].Points.AddXY(x[i], y[i]);
            }
            for (int i = 0; i < nbig; i++)//вывод первой функции 
            {
                chart1.Series[1].Points.AddXY(x[i], s[i]);
            }

        }
        private void massiv2()
        {
            double h = Convert.ToDouble(textBox1.Text);
            double hsmall = Convert.ToDouble(textBox2.Text);
            int n = Convert.ToInt32(2 / h) + 1;
            int nbig = Convert.ToInt32(2 / hsmall) + 1;
            double[] x = new double[nbig];
            double[] y = new double[nbig];
            double[] a = new double[n];
            double[] b = new double[n];
            double[] c = new double[n];
            double[] d = new double[n];
            double[] u = new double[n];
            double[] v = new double[n];
            double[] m = new double[n];
            double[] f = new double[n];
            double[] xi = new double[n];
            double[] s = new double[nbig+200];
            xi[0] = 0;
            for (int i = 1; i < n; i++)
            {
                xi[i] = xi[i - 1] + h;
            }

            for (int i = 0; i < n; i++)
            {
                f[i] = Math.Exp(-xi[i])*xi[i];//узлы
            }

            x[0] = 0;
            for (int i = 1; i < nbig; i++)
            {
                x[i] = x[i - 1] + hsmall;
            }

            for (int i = 0; i < nbig; i++)//сама функция
            {
                y[i] = Math.Exp(-x[i]) * x[i];
            }

            a[0] = 1;
            a[n - 1] = 1;
            b[0] = 0;
            c[0] = 0;
            c[n - 1] = 0;
            d[0] = -2 * Math.Exp(-xi[0]) + x[0] * Math.Exp(-xi[0]);
            d[n - 1] = -2 * Math.Exp(-xi[n-1]) + x[n-1] * Math.Exp(-xi[n-1]);
            for (int i = 1; i < n - 1; i++)//матрица для m 
            {
                a[i] = 4;
                b[i] = 1;
                c[i] = 1;
                d[i] = 6 * (f[i + 1] - 2 * f[i] + f[i - 1]) / (h * h);
            }

            u[0] = b[0] / a[0] * (-1);
            v[0] = d[0] / a[0];
            for (int i = 1; i < n; i++)//вспомогательные переменные для m 
            {
                u[i] = (-1) * b[i] / (c[i] * u[i - 1] + a[i]);
                v[i] = (d[i] - c[i] * v[i - 1]) / (c[i] * u[i - 1] + a[i]);
            }

            m[n - 1] = v[n - 1];
            for (int i = n - 2; i > -1; i--)//вычисление m 
            {
                m[i] = u[i] * m[i + 1] + v[i];
            }
            int k = 0;
            double t;
            for (int i = 0; i < n - 1; i++)
            {
                for (int j = 0; j < h / hsmall; j++)
                {
                    t = hsmall * j / h;
                    s[k] = f[i] * psi1(t) + f[i + 1] * t + m[i] * h * h * psi3(t) + m[i + 1] * h * h * psi4(t);
                    k++;
                }
            }
            t = 1;
            s[nbig - 1] = f[n - 2] * psi1(t) + f[n - 1] * t + m[n - 2] * h * h * psi3(t) + m[n - 1] * h * h * psi4(t);

            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            for (int i = 0; i < nbig; i++)//вывод первой функции 
            {
                chart1.Series[0].Points.AddXY(x[i], y[i]);
            }
            for (int i = 0; i < nbig; i++)//вывод первой функции 
            {
                chart1.Series[1].Points.AddXY(x[i], s[i]);
            }

        }
        private void massiv3()
        {
            double h = Convert.ToDouble(textBox1.Text);
            double hsmall = Convert.ToDouble(textBox2.Text);
            int n = Convert.ToInt32(2*Math.PI / h) + 1;
            int nbig = Convert.ToInt32(2 * Math.PI / hsmall) + 1;
            double[] x = new double[nbig];
            double[] y = new double[nbig];
            double[] a = new double[n];
            double[] b = new double[n];
            double[] c = new double[n];
            double[] d = new double[n];
            double[] u = new double[n];
            double[] v = new double[n];
            double[] m = new double[n];
            double[] f = new double[n];
            double[] xi = new double[n];
            double[] s = new double[nbig+1000];
            xi[0] = -Math.PI;
            for (int i = 1; i < n; i++)
            {
                xi[i] = xi[i - 1] + h;
            }

            for (int i = 0; i < n; i++)
            {
                f[i] = Math.Sin(xi[i]);//узлы
            }

            x[0] = -Math.PI;
            for (int i = 1; i < nbig; i++)
            {
                x[i] = x[i - 1] + hsmall;
            }

            for (int i = 0; i < nbig; i++)//сама функция
            {
                y[i] = Math.Sin(x[i]);
            }

            a[0] = 1;
            a[n - 1] = 1;
            b[0] = 0;
            c[0] = 0;
            c[n - 1] = 0;
            d[0] = -Math.Sin(xi[0]);
            d[n - 1] = -Math.Sin(xi[n-1]);
            for (int i = 1; i < n - 1; i++)//матрица для m 
            {
                a[i] = 4;
                b[i] = 1;
                c[i] = 1;
                d[i] = 6 * (f[i + 1] - 2 * f[i] + f[i - 1]) / (h * h);
            }

            u[0] = b[0] / a[0] * (-1);
            v[0] = d[0] / a[0];
            for (int i = 1; i < n; i++)//вспомогательные переменные для m 
            {
                u[i] = (-1) * b[i] / (c[i] * u[i - 1] + a[i]);
                v[i] = (d[i] - c[i] * v[i - 1]) / (c[i] * u[i - 1] + a[i]);
            }

            m[n - 1] = v[n - 1];
            for (int i = n - 2; i > -1; i--)//вычисление m 
            {
                m[i] = u[i] * m[i + 1] + v[i];
            }
            int k = 0;
            double t;
            for (int i = 0; i < n - 1; i++)
            {
                for (int j = 0; j < h / hsmall; j++)
                {
                    t = hsmall * j / h;
                    s[k] = f[i] * psi1(t) + f[i + 1] * t + m[i] * h * h * psi3(t) + m[i + 1] * h * h * psi4(t);
                    k++;
                }
            }
            

            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            for (int i = 0; i < nbig; i++)//вывод первой функции 
            {
                chart1.Series[0].Points.AddXY(x[i], y[i]);
            }
            for (int i = 0; i < nbig; i++)//вывод первой функции 
            {
                chart1.Series[1].Points.AddXY(x[i], s[i]);
            }

        }
        private double psi1(double t)
        {
            t = 1 - t;
            return t;
        }

        private double psi3(double t)
        {
            t = -1 * t * (t - 1) * (t - 2) / 6;
            return t;
        }

        private double psi4(double t)
        {
            t = t * (t - 1) * (t + 1) / 6;
            return t;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            //massiv1();
            
            //шаг большой 
            double h = Convert.ToDouble(textBox1.Text);
            int N = Convert.ToInt32((2 / h) + 1);
            double[] x = new double[N];

            //шаг маленький 
            double hm = Convert.ToDouble(textBox2.Text);
            int Nm = Convert.ToInt32((2 / hm) + 1);

            //график по ответу 
            //задаём x и y 
            double[] xm = new double[Nm];
            double[] y1 = new double[Nm];
            for (int i = 0; i < Nm; i++)
            {
                xm[i] = -1 + Convert.ToDouble(i) * hm;
                y1[i] = Math.Exp(-Math.Pow(xm[i], 2));
            }

            //функция в узлах: 
            double[] y = new double[N];
            for (int i = 0; i < N; i++)
            {
                y[i] = Math.Exp(-Math.Pow(-1 + Convert.ToDouble(i) * h, 2));
            }

            // ПРОГОНКА 
           
            //массивы a,b,c,d 
            double[] a = new double[N];
            a[0] = 1; a[N - 1] = 1;
            for (int i = 1; i < a.Length - 1; i++)
            {
                a[i] = 4;
            }

            double[] b = new double[N];
            b[0] = 0; b[N - 1] = 0;
            for (int i = 1; i < b.Length - 1; i++)
            {
                b[i] = 1;
            }

            double[] c = new double[N];
            c[0] = 0; 
            for (int i = 1; i < c.Length - 1; i++)
            {
                c[i] = 1;
            }
            c[N - 1] = 0;
            double[] d = new double[N];
            d[0] = -2 / Math.Exp(xm[0] * xm[0]) + 4 * xm[0] * xm[0] / Math.Exp(xm[0] * xm[0]);
            d[N - 1] = -2 / Math.Exp(xm[Nm - 1] * xm[Nm - 1]) + 4 * xm[Nm - 1] * xm[Nm - 1] / Math.Exp(xm[Nm - 1] * xm[Nm - 1]);
            for (int i = 1; i < d.Length - 1; i++)
            {
                d[i] = 6 * (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
            }
 
            //прямой ход 
            //задаём U и V 
            double[] U = new double[N];
            U[0] = b[0] / a[0]*(-1); //U[N - 1] = 0; 
            for (int j = 1; j < N - 1; j++)
            {
                U[j] = -b[j] / (c[j] * U[j - 1] + a[j]);
            }
            double[] V = new double[N];
            V[0] = d[0] / a[0];
            for (int j = 1; j < N; j++)
            {
                V[j] = (d[j] - c[j] * V[j - 1]) / (c[j] * U[j - 1] + a[j]);
            }

            //обратный ход 
            double[] m = new double[N];
            m[N - 1] = V[N - 1];
            for (int i = N - 2; i >= 0; i--)
            {
                m[i] = U[i] * m[i + 1] + V[i];
            }
           
            //double j= y[i - 1] + y[i + 1] - y[] * 2
            //СПЛАЙН 
           

            double t; double[] S = new double[Nm + 300];

            int iks = 0;
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = 0; j < h / hm; j++)
                {
                    t = j * hm / h;
                    S[iks] = y[i] * Fi1(t) + y[i + 1] * t + m[i] * h * h * Fi3(t) + m[i + 1] * h * h * Fi4(t);
                    iks++;
                }
            }
            //вычисление в последней точке 
            if (iks < Nm + 1)
            {
                t = 1;
                S[Nm - 1] = y[N - 2] * Fi1(t) + y[N - 1] * t + m[N - 2] * h * h * Fi3(t) + m[N - 1] * h * h * Fi4(t);
            }
            //рисуем 
            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            for (int i = 0; i < Nm; i++)//по ответу 
                chart1.Series[0].Points.AddXY(xm[i], y1[i]);
            for (int i = 0; i < Nm; i++)//по сплайнам 
                chart1.Series[1].Points.AddXY(xm[i], S[i]);

        }



        //Fi 
        private double Fi1(double t)
        {
            double f1 = 1 - t;
            return f1;
        }
        private double Fi3(double t)
        {
            double f3 = (-1) * t * (t - 1) * (t - 2) / 6;
            return f3;
        }
        private double Fi4(double t)
        {
            double f4 = 1 * t * (t - 1) * (t + 1) / 6;
            return f4;
        }

        private void button2_Click(object sender, EventArgs e) // 2 - для функции exp * x 
        {
            //шаг большой 
            double h = Convert.ToDouble(textBox1.Text);
            int N = Convert.ToInt32((2 / h) + 1);
            double[] x = new double[N];

            //шаг маленький 
            double hm = Convert.ToDouble(textBox2.Text);
            int Nm = Convert.ToInt32((2 / hm) + 1);

            //график по ответу 
            //задаём x и y 
            double[] xm = new double[Nm];
            double[] y1 = new double[Nm];
            for (int i = 0; i < Nm; i++)
            {
                xm[i] = Convert.ToDouble(i) * hm;
                y1[i] = Math.Exp(-xm[i]) * xm[i];
            }

            //функция в узлах: 
            double[] y = new double[N];
            for (int i = 0; i < N; i++)
            {
                y[i] = Math.Exp(-Convert.ToDouble(i) * h) * Convert.ToDouble(i) * h;
            }

            // ПРОГОНКА 

            //массивы a,b,c,d 
            double[] a = new double[N];
            a[0] = 1; a[N - 1] = 1;
            for (int i = 1; i < a.Length - 1; i++)
            {
                a[i] = 4;
            }

            double[] b = new double[N];
            b[0] = 0; b[N - 1] = 0;
            for (int i = 1; i < b.Length - 1; i++)
            {
                b[i] = 1;
            }

            double[] c = new double[N];
            c[1] = 1; c[N - 1] = 0;
            for (int i = 2; i < c.Length - 1; i++)
            {
                c[i] = 1;
            }

            double[] d = new double[N];
            d[0] = Math.Exp(-xm[0]) * (xm[0] - 2);
            d[N - 1] = Math.Exp(-xm[Nm - 1]) * (xm[Nm - 1] - 2);
            for (int i = 1; i < d.Length - 1; i++)
            {
                d[i] = 6 * (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
            }

            //прямой ход 
            //задаём U и V 
            double[] U = new double[N];
            U[0] = -b[0] / a[0]; U[N - 1] = 0;
            for (int j = 1; j < N - 1; j++)
            {
                U[j] = -b[j] / (c[j] * U[j - 1] + a[j]);
            }
            double[] V = new double[N];
            V[0] = d[0] / a[0];
            for (int j = 1; j < N; j++)
            {
                V[j] = (d[j] - c[j] * V[j - 1]) / (c[j] * U[j - 1] + a[j]);
            }

            //обратный ход 
            double[] m = new double[N];
            m[N - 1] = V[N - 1];
            for (int i = N - 2; i >= 0; i--)
            {
                m[i] = U[i] * m[i + 1] + V[i];
            }

            //СПЛАЙН 

            double t; double[] S = new double[Nm + 300];
            int iks = 0;
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = 0; j < h / hm; j++)
                {
                    t = j * hm / h;
                    S[iks] = y[i] * Fi1(t) + y[i + 1] * t + m[i] * h * h * Fi3(t) + m[i + 1] * h * h * Fi4(t);
                    iks++;
                }
            }
            //вычисление в последней точке 
            if (iks < Nm + 1)
            {
                t = 1;
                S[Nm - 1] = y[N - 2] * Fi1(t) + y[N - 1] * t + m[N - 2] * h * Fi3(t) + m[N - 1] * h * Fi4(t);
            }
            //рисуем 
            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            for (int i = 0; i < Nm; i++)//по ответу 
                chart1.Series[0].Points.AddXY(xm[i], y1[i]);
            for (int i = 0; i < Nm; i++)//по сплайнам 
                chart1.Series[1].Points.AddXY(xm[i], S[i]);

        }



        private void button3_Click(object sender, EventArgs e) // 3 - для функции sin 
        {
            //шаг большой 
            double h = Convert.ToDouble(textBox1.Text);
            int N = Convert.ToInt32((2 * Math.PI / h) + 1);
            double[] x = new double[N];

            //шаг маленький 
            double hm = Convert.ToDouble(textBox2.Text);
            int Nm = Convert.ToInt32((2 * Math.PI / hm) + 1);

            //график по ответу 
            //задаём x и y 
            double[] xm = new double[Nm];
            double[] y1 = new double[Nm];
            for (int i = 0; i < Nm; i++)
            {
                xm[i] = -Math.PI + Convert.ToDouble(i) * hm;
                y1[i] = Math.Sin(xm[i]);
            }

            //функция в узлах: 
            double[] y = new double[N];
            for (int i = 0; i < N; i++)
            {
                y[i] = Math.Sin(-Math.PI + Convert.ToDouble(i) * h);
            }

            // ПРОГОНКА 

            //массивы a,b,c,d 
            double[] a = new double[N];
            a[0] = 1; a[N - 1] = 1;
            for (int i = 1; i < a.Length - 1; i++)
            {
                a[i] = 4;
            }

            double[] b = new double[N];
            b[0] = 0; b[N - 1] = 0;
            for (int i = 1; i < b.Length - 1; i++)
            {
                b[i] = 1;
            }

            double[] c = new double[N];
            c[1] = 1; c[N - 1] = 0;
            for (int i = 2; i < c.Length - 1; i++)
            {
                c[i] = 1;
            }

            double[] d = new double[N];
            d[0] = -Math.Sin(xm[0]);
            d[N - 1] = -Math.Sin(xm[Nm - 1]);
            for (int i = 1; i < d.Length - 1; i++)
            {
                d[i] = 6 * (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
            }

            //прямой ход 
            //задаём U и V 
            double[] U = new double[N];
            U[0] = -b[0] / a[0]; U[N - 1] = 0;
            for (int j = 1; j < N - 1; j++)
            {
                U[j] = -b[j] / (c[j] * U[j - 1] + a[j]);
            }
            double[] V = new double[N];
            V[0] = d[0] / a[0];
            for (int j = 1; j < N; j++)
            {
                V[j] = (d[j] - c[j] * V[j - 1]) / (c[j] * U[j - 1] + a[j]);
            }

            //обратный ход 
            double[] m = new double[N];
            m[N - 1] = V[N - 1];
            for (int i = N - 2; i >= 0; i--)
            {
                m[i] = U[i] * m[i + 1] + V[i];
            }

            //СПЛАЙН 

            double t; double[] S = new double[Nm + 1000];
            int iks = 0;
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = 0; j < h / hm; j++)
                {
                    t = j * hm / h;
                    S[iks] = y[i] * Fi1(t) + y[i + 1] * t + m[i] * h * h * Fi3(t) + m[i + 1] * h * h * Fi4(t);
                    iks++;
                }
            }
           double[] prov1 = new double[N];
           double[] prov2 = new double[N];
           for (int i = 1; i < N - 1; i++)
           {
               prov1[i] = m[i - 1] + 4 * m[i] + m[i + 1];
               prov2[i] = 6 * (y[i - 1] + y[i + 1] - y[i] * 2) / (h * h);
           }
            //вычисление в последней точке 
            if (iks < Nm + 1)
            {
                t = 1;
                S[Nm - 1] = y[N - 2] * Fi1(t) + y[N - 1] * t + m[N - 2] * h * h * Fi3(t) + m[N - 1] * h * h * Fi4(t);
            }
            //рисуем 
            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            for (int i = 0; i < Nm; i++)//по ответу 
                chart1.Series[0].Points.AddXY(xm[i], y1[i]);
            for (int i = 0; i < Nm; i++)//по сплайнам 
                chart1.Series[1].Points.AddXY(xm[i], S[i]);
        }
    }

}

/*private void button2_Click(object sender, EventArgs e)
        {
            massiv2();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            massiv3();
        }
    }
}*/
