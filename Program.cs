using System;
using System.Collections.Generic;
//using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;

using System.Diagnostics; // подключаем пространство имён, в котором находится класс Stopwatch
using System.Threading;

namespace ConsoleApplication1
{

    class Program
    {
        static public double[] k(double t, double[] x, double h)
        {
            double[] x1 = new double[6];
            double[] k1 = new double[6];
            double[] k2 = new double[6];
            double[] k3 = new double[6];
            double[] k4 = new double[6];
            double[] f = new double[6];
            double[] xx = new double[6];
            //1
            f = fcn(t, x);
            k1 = f;
            //2
            for (int i = 0; i < 6; i++) { xx[i] = x[i] + h * k1[i] / 2; }
            f = fcn(t + h / 2, xx);
            k2 = f;
            //3
            for (int i = 0; i < 6; i++) { xx[i] = x[i] + h * k2[i] / 2; }
            f = fcn(t + h / 2, xx);
            k3 = f;
            //4
            for (int i = 0; i < 6; i++) { xx[i] = x[i] + h * k3[i]; }
            f = fcn(t + h, xx);
            k4 = f;
            //пересчет
            for (int i = 0; i < 6; i++) { x1[i] = x[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6; }
            return x1;
        }

        static public double[] fcn(double t, double[] x)
        {
            double[] f = new double[6];
            double r;
            double mu = 2.959122082855911025E-4;
            r = Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            f[0] = x[3];
            f[1] = x[4];
            f[2] = x[5];
            f[3] = -mu / r / r / r * x[0];
            f[4] = -mu / r / r / r * x[1];
            f[5] = -mu / r / r / r * x[2];
            return f;
        }

        static public double dr(double[] x0, double[] x)
        {
            double d;
            d = Math.Sqrt((x[0] - x0[0]) * (x[0] - x0[0]) + (x[1] - x0[1]) * (x[1] - x0[1]) + (x[2] - x0[2]) * (x[2] - x0[2]));
            return d;
        }



        static void Main(string[] args)
        {
            Stopwatch time = new Stopwatch(); // создаём объект Stopwatch
            time.Start(); // запускаем отсчёт времени
            //создаю файлы для записи данных
            System.IO.StreamWriter t1 = new System.IO.StreamWriter(@"C:\\t.txt");
            System.IO.StreamWriter t2 = new System.IO.StreamWriter(@"C:\\H.txt");
            System.IO.StreamWriter t3 = new System.IO.StreamWriter(@"C:\\R.txt");
            System.IO.StreamWriter t4 = new System.IO.StreamWriter(@"C:\\Otvet.txt");
            double hpred;
            //-----------------------------------//
            double[] x = new double[6];
            double[] x0 = new double[6];
            double[] x1 = new double[6];
            double[] x2 = new double[6];
            double ts, tf, h, t, e_cal, d, r;
            double e_tol = 1E-8;
            //!body
            ts = 2457700.5;
            tf = 2459380.5;
            x[0] = 2.42018325897182938;  //cer
            x[1] = 1.52163771807878472;  //cer
            x[2] = 0.22441527622794077;  //cer
            x[3] = -0.00558516807878721; //cer
            x[4] = 0.00697766230356442;  //cer
            x[5] = 0.00442730669383955;  //cer

            x0[0] = 2.42018325897182938;  //cer
            x0[1] = 1.52163771807878472;  //cer
            x0[2] = 0.22441527622794077;  //cer
            x0[3] = -0.00558516807878721; //cer
            x0[4] = 0.00697766230356442;  //cer
            x0[5] = 0.00442730669383955;  //cer

            //tf = 2458223.50;            //phaeton 
            //x[0] = 1.2955150372595632;  //phaeton
            //x[1] = 0.43730519196972016; //phaeton
            //x[2] = 0.73949230027487296; //phaeton
            //x[3] = 0.0066601612740838;  //phaeton
            //x[4] = 0.00813247761290623; //phaeton
            //x[5] = 0.006119327183538;   //phaeton

            //x0[0] = 1.2955150372595632;  //phaeton
            //x0[1] = 0.43730519196972016; //phaeton
            //x0[2] = 0.73949230027487296; //phaeton
            //x0[3] = 0.0066601612740838;  //phaeton
            //x0[4] = 0.00813247761290623; //phaeton
            //x0[5] = 0.006119327183538;   //phaeton
            h = 0.001;
            t = ts;
            while (t < tf)
            {
                hpred = h;
                x1 = k(t, x, h);
                x2 = k(t, x, h / 2);
                x2 = k(t + h / 2, x2, h / 2);
                d = dr(x1, x2);
                e_cal = d / Math.Pow((1 - 0.5), 0.2);//posmotret formyly
                t = t + h;
                x = x2;
                h = h * Math.Pow((e_tol / e_cal), 0.2);//posmotret formyly
                r = Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
                if (d == 0) { h = hpred; }
                //write t,h,r
                t1.WriteLine(t);
                t2.WriteLine(h);
                t3.WriteLine(r);
                System.Console.WriteLine(t);
                if ((((t + h) > tf)) & (Math.Abs(t - tf) > e_tol)) { h = tf - t; }
            }
            //teper schitaem velichiny oshibki
            //write t,x
            t4.WriteLine(t);
            t4.WriteLine(x[0]); t4.WriteLine(x[1]); t4.WriteLine(x[2]); t4.WriteLine(x[3]); t4.WriteLine(x[4]); t4.WriteLine(x[5]);
            h = -h;
            while (t > ts)
            {
                hpred = h;
                x1 = k(t, x, h);
                x2 = k(t, x, h / 2);
                x2 = k(t + h / 2, x2, h / 2);
                d = dr(x1, x2);
                e_cal = d / Math.Pow((1 - 0.5), 0.2);//posmotret formyly
                t = t + h;
                x = x2;
                h = h * Math.Pow((e_tol / e_cal), 0.2);//posmotret formyly
                if (d == 0) { h = hpred; }
                if (((t + h) < ts) & (Math.Abs(t - ts) > e_tol)) { h = ts - t; }
            }
            d = dr(x, x0);
            //write dr=d
            t4.WriteLine(d);
            t4.WriteLine(time.Elapsed); // выводим затраченное время
            //end program Runge_Kutt
            t1.Close();
            t2.Close();
            t3.Close();
            t4.Close();
            Console.ReadLine();
        }
    }
}