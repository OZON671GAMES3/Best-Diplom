using System;
using System.Collections.Generic;
//using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;

namespace ConsoleApplication1
{
    //ОПИСАНИЕ КЛАССА ДРЕДЕР
    class DEreader
    {

        class DEconst
        {
            public double AE;
            public double speedLight;
            public double G;
            public double rSun;
            public double compSun;
            public double rEarth;
            public double compEarth;
            public double[] massPlanet;
            public double massMoon;
            public double[] massAst;


            public DEconst()
            {
                massPlanet = new double[9];//[9]
                massAst = new double[3];//[3]
            }
        };

        DEconst deConst = new DEconst();

        public
            class DEheader
        {
            public int sizeStr;
            public int step;
            public double tMin;
            public double tMax;
            public double EarthDivMoon;
            public int[] nper = new int[11];//[11]
            public int[] pow = new int[11];//[11]
            public double[] raz = new double[11];//[11]
        };

        DEheader deHeader = new DEheader();


        List<double> buf = new List<double>();
        String sPathFond;
        double dJD;

        //КОНСТРУКТОРЫ КЛАССА ДРЭДЕР
        public DEreader()
        {
        }

        /*----------------------------------------------------------------------------*/
        public DEreader(int fond, String pathFond)
        {
            sPathFond = pathFond;
            if (fond == 405)
            {
                deHeader.sizeStr = 1018;
                deHeader.step = 32;
                deHeader.tMin = 2305424.50;
                deHeader.tMax = 2525008.50;
                deHeader.EarthDivMoon = 0.813005600000000044E+02;

                deHeader.nper[0] = 3;
                deHeader.nper[1] = 171;
                deHeader.nper[2] = 231;
                deHeader.nper[3] = 309;
                deHeader.nper[4] = 342;
                deHeader.nper[5] = 366;
                deHeader.nper[6] = 387;
                deHeader.nper[7] = 405;
                deHeader.nper[8] = 423;
                deHeader.nper[9] = 441;
                deHeader.nper[10] = 753;

                deHeader.pow[0] = 14;
                deHeader.pow[1] = 10;
                deHeader.pow[2] = 13;
                deHeader.pow[3] = 11;
                deHeader.pow[4] = 8;
                deHeader.pow[5] = 7;
                deHeader.pow[6] = 6;
                deHeader.pow[7] = 6;
                deHeader.pow[8] = 6;
                deHeader.pow[9] = 13;
                deHeader.pow[10] = 11;

                deHeader.raz[0] = deHeader.step / 4;
                deHeader.raz[1] = deHeader.step / 2;
                deHeader.raz[2] = deHeader.step / 2;
                deHeader.raz[3] = deHeader.step / 1;
                deHeader.raz[4] = deHeader.step / 1;
                deHeader.raz[5] = deHeader.step / 1;
                deHeader.raz[6] = deHeader.step / 1;
                deHeader.raz[7] = deHeader.step / 1;
                deHeader.raz[8] = deHeader.step / 1;
                deHeader.raz[9] = deHeader.step / 8;
                deHeader.raz[10] = deHeader.step / 2;

                deConst.AE = 149597870.691000015;
                deConst.speedLight = 299792.457999999984;
                deConst.G = 0.295912208285591095E-03;
                deConst.rSun = 696000.000000000000;
                deConst.compSun = 0.199999999999999991E-06;
                deConst.rEarth = 6378.13699999999972;
                deConst.compEarth = 0.108262599999999994E-02;

                deConst.massPlanet[0] = 0.491254745145081187E-10;
                deConst.massPlanet[1] = 0.724345248616270270E-09;
                deConst.massPlanet[2] = 0.899701134671249882E-09 * (1 - 1 / (deHeader.EarthDivMoon + 1)); // Earth
                deConst.massPlanet[3] = 0.954953510577925806E-10;
                deConst.massPlanet[4] = 0.282534590952422643E-06;
                deConst.massPlanet[5] = 0.845971518568065874E-07;
                deConst.massPlanet[6] = 0.129202491678196939E-07;
                deConst.massPlanet[7] = 0.152435890078427628E-07;
                deConst.massPlanet[8] = 0.218869976542596968E-11;

                deConst.massMoon = 0.899701134671249882E-09 / (deHeader.EarthDivMoon + 1);

                deConst.massAst[0] = 0.139078737894227800E-12;
                deConst.massAst[1] = 0.295912208285591104E-13;
                deConst.massAst[2] = 0.384685870771268372E-13;
            }
            dJD = 0;
        }

        ~DEreader()
        {
            buf.Clear();
        }

        /*----------------------------------------------------------------------------*/
        public void DEreaderGetPlanetPoz(bool vel, double jdate, int index, bool geleo, double[] poz)
        {

            if (jdate != dJD)
            {
                read(jdate);
                dJD = jdate;
            }
            int max;
            if (vel) { max = 6; } else max = 3;
            double[] tc = new double[16];//[16]
            double[] tcp = new double[16];//[16]
            for (int i = 0; i < 16; i++)
            {
                tc[i] = 0;
                tcp[i] = 0;
            }

            int n1 = 0;
            double d = buf[0];
            double c = d - deHeader.raz[index];
            n1 = -3 * deHeader.pow[index];
            do
            {
                c += deHeader.raz[index];
                d += deHeader.raz[index];
                n1 += 3 * deHeader.pow[index];
            } while (jdate < c || jdate > d);
            cheb(vel, c, d, jdate, deHeader.pow[index], tc, tcp);
            coor(vel, index, n1, tc, tcp, poz);

            if (index == 2 || index == 9)
            {
                int idx = index == 2 ? 9 : 2;
                n1 = 0;
                d = buf[0];
                c = d - deHeader.raz[idx];
                n1 = -3 * deHeader.pow[idx];
                do
                {
                    c += deHeader.raz[idx];
                    d += deHeader.raz[idx];
                    n1 += 3 * deHeader.pow[idx];
                } while (jdate < c || jdate > d);
                cheb(vel, c, d, jdate, deHeader.pow[idx], tc, tcp);
                double[] x = new double[max];//[6]
                coor(vel, idx, n1, tc, tcp, x);
                if (index == 2)
                    for (int k = 0; k < max; k++)
                        poz[k] = poz[k] - x[k] / (deHeader.EarthDivMoon + 1); // земля  
                if (index == 9)
                    for (int k = 0; k < max; k++)
                        x[k] = x[k] - poz[k] / (deHeader.EarthDivMoon + 1); // земля
                if (index == 9)
                    for (int k = 0; k < max; k++)
                        poz[k] += x[k]; // луна
            }
            if (geleo)
            {
                n1 = 0;
                d = buf[0];
                c = d - deHeader.raz[10];
                n1 = -3 * deHeader.pow[10];
                do
                {
                    c += deHeader.raz[10];
                    d += deHeader.raz[10];
                    n1 += 3 * deHeader.pow[10];
                } while (jdate < c || jdate > d);
                cheb(vel, c, d, jdate, deHeader.pow[10], tc, tcp);
                double[] x = new double[max];//[6]
                coor(vel, 10, n1, tc, tcp, x);

                for (int k = 0; k < max; k++)
                    poz[k] = (poz[k] - x[k]) / deConst.AE;
            }
            else
                for (int k = 0; k < max; k++)
                    poz[k] /= deConst.AE;
        }

        void cheb(bool vel, double a, double b, double t, int st, double[] tc, double[] tcp)
        {
            double rat;
            double dlin = b - a;
            double tau = 2 * (t - a) / dlin - 1;
            double tau2 = tau * 2;
            tc[0] = 1;
            tc[1] = tau;
            for (int i = 2; i <= st; i++)
                tc[i] = tau2 * tc[i - 1] - tc[i - 2];
            if (vel)
            {
                rat = 4 / dlin;
                tcp[0] = 0;
                tcp[1] = 2 / dlin;
                for (int i = 2; i <= st; i++)
                    tcp[i] = rat * tc[i - 1] + tau2 * tcp[i - 1] - tcp[i - 2];
            }
        }

        /*----------------------------------------------------------------------------*/
        void coor(bool vel, int i, int n1, double[] tc, double[] tcp, double[] x)
        {
            int jj;
            double s;
            int ist = deHeader.pow[i];
            for (int k = 0; k < 3; k++)
            {
                s = 0;
                jj = deHeader.nper[i] + k * ist + n1 - 1;
                for (int j = 0; j < deHeader.pow[i]; j++)
                    s += tc[j] * buf[jj + j];
                x[k] = s;
            }
            if (vel)
                for (int k = 0; k < 3; k++)
                {
                    s = 0;
                    jj = deHeader.nper[i] + k * ist + n1 - 1;
                    s = 0;
                    for (int j = 0; j < deHeader.pow[i]; j++)
                        s += tcp[j] * buf[jj + j];
                    x[k + 3] = s;
                }
        }

        /*----------------------------------------------------------------------------*/
        void read(double t)
        {
            int nrc = (int)Math.Truncate((dJD - deHeader.tMin) / deHeader.step);
            int nr = (int)Math.Truncate((t - deHeader.tMin) / deHeader.step);

            if (nrc != nr)
            {
                buf.Clear();
                var fond = File.Open(sPathFond, FileMode.Open);
                var binfond = new BinaryReader(fond);
                fond.Seek(nr * deHeader.sizeStr * sizeof(double), SeekOrigin.Begin);
                for (int i = 0; i < deHeader.sizeStr; i++)
                {
                    double d;
                    d = binfond.ReadDouble();
                    buf.Add(d);
                }
                fond.Close();

            }
        }
    }


    //-----------------------------------------------------------------------------------Старое
    class Program
    {
        public class Vector3
        {
            public double[] values = new double[3];

            public double _1
            {

                get { return values[0]; }
                set { values[0] = value; }
            }
            public double _2
            {
                get { return values[1]; }
                set { values[1] = value; }
            }
            public double _3
            {
                get { return values[2]; }
                set { values[2] = value; }
            }

            public Vector3 Clone()
            {
                Vector3 v = new Vector3();
                v.values = (double[])values.Clone();
                return v;
            }
        }

        class Coordinates
        {
            public Vector3 pos = new Vector3(), v = new Vector3();
            public Coordinates Clone()
            {
                Coordinates c = new Coordinates();
                c.pos = pos.Clone();
                c.v = v.Clone();
                return c;
            }

        }

        const double G = 2.959122082855911025E-4;
        static double H = 0.001;
        const double EpsPresset = 1E-14;


        static string filePath = "C:\\Users\\Amanda\\Desktop\\coord.txt";
        static string output = "C:\\Users\\Amanda\\Desktop\\newcoord.txt";
        static string dora = "C:\\Users\\Amanda\\Desktop\\H.txt";
        /////  что то с ними не то
        public static double Koordinata(double V)
        {
            return V;
        }

        public static double Skorost(double x, double R)
        {
            return (-G * x) / (R * R * R);
        }

        static double calcR(Coordinates input, Coordinates k, double h)
        {
            return Math.Sqrt(Math.Pow((input.pos._1 + k.pos._1 * h), 2) + Math.Pow((input.pos._2 + k.pos._2 * h), 2) + Math.Pow((input.pos._3 + k.pos._3 * h), 2));
        }
        // Основной!!!
        static Coordinates calc(Coordinates input, double t0, double tn, DEreader DE405, double jdate)
        {
            double t = t0;
            double halfH = H / 2;
            Coordinates output = input.Clone();
            Coordinates k = new Coordinates(), k2 = new Coordinates(), k3 = new Coordinates(), k4 = new Coordinates();
            Coordinates H2;
            double qqq = 0;
            double R = Math.Sqrt(Math.Pow(output.pos._1, 2) + Math.Pow(output.pos._2, 2) + Math.Pow(output.pos._3, 2));
            ///// начало цикла
            System.IO.StreamWriter textFile = new System.IO.StreamWriter(@"C:\\shag.txt");
            System.IO.StreamWriter File = new System.IO.StreamWriter(@"C:\\radius.txt");
            System.IO.StreamWriter text = new System.IO.StreamWriter(@"C:\\t.txt");
            while (t < tn)
            {
                if (t + H > tn)
                { H = tn - t; };

                halfH = H / 2.0;
                k.pos._1 = Koordinata(output.v._1);   ///Координаты
                k.pos._2 = Koordinata(output.v._2);
                k.pos._3 = Koordinata(output.v._3);

                k.v._1 = Skorost(output.pos._1, R) + Vozmysheniya(t, 0, output.pos, DE405, jdate);       /// Скорости
                k.v._2 = Skorost(output.pos._2, R) + Vozmysheniya(t, 1, output.pos, DE405, jdate);
                k.v._3 = Skorost(output.pos._3, R) + Vozmysheniya(t, 2, output.pos, DE405, jdate);

                R = calcR(output, k, halfH);
                //////
                k2.pos._1 = Koordinata(output.v._1 + k.v._1 * halfH);
                k2.pos._2 = Koordinata(output.v._2 + k.v._2 * halfH);
                k2.pos._3 = Koordinata(output.v._3 + k.v._3 * halfH);

                k2.v._1 = Skorost(output.pos._1 + k.pos._1 * halfH, R) + Vozmysheniya(t, 0, k.pos, DE405, jdate);
                k2.v._2 = Skorost(output.pos._2 + k.pos._2 * halfH, R) + Vozmysheniya(t, 1, k.pos, DE405, jdate);
                k2.v._3 = Skorost(output.pos._3 + k.pos._3 * halfH, R) + Vozmysheniya(t, 2, k.pos, DE405, jdate);

                R = calcR(output, k2, halfH);
                ///////
                k3.pos._1 = Koordinata(output.v._1 + k2.v._1 * halfH);
                k3.pos._2 = Koordinata(output.v._2 + k2.v._2 * halfH);
                k3.pos._3 = Koordinata(output.v._3 + k2.v._3 * halfH);

                k3.v._1 = Skorost(output.pos._1 + k2.pos._1 * halfH, R) + Vozmysheniya(t, 0, k2.pos, DE405, jdate);
                k3.v._2 = Skorost(output.pos._2 + k2.pos._2 * halfH, R) + Vozmysheniya(t, 1, k2.pos, DE405, jdate);
                k3.v._3 = Skorost(output.pos._3 + k2.pos._3 * halfH, R) + Vozmysheniya(t, 2, k2.pos, DE405, jdate);

                R = calcR(output, k3, H);
                /////////
                k4.pos._1 = Koordinata(output.v._1 + k3.v._1 * H);
                k4.pos._2 = Koordinata(output.v._2 + k3.v._2 * H);
                k4.pos._3 = Koordinata(output.v._3 + k3.v._3 * H);

                k4.v._1 = Skorost(output.pos._1 + k3.pos._1 * H, R) + Vozmysheniya(t, 0, k3.pos, DE405, jdate);
                k4.v._2 = Skorost(output.pos._2 + k3.pos._2 * H, R) + Vozmysheniya(t, 1, k3.pos, DE405, jdate);
                k4.v._3 = Skorost(output.pos._3 + k3.pos._3 * H, R) + Vozmysheniya(t, 2, k3.pos, DE405, jdate);

                /// Перещет скоростей и координат
                output.pos._1 = recalculation(output.pos._1, k.pos._1, k2.pos._1, k3.pos._1, k4.pos._1, H / 6.0);
                output.pos._2 = recalculation(output.pos._2, k.pos._2, k2.pos._2, k3.pos._2, k4.pos._2, H / 6.0);
                output.pos._3 = recalculation(output.pos._3, k.pos._3, k2.pos._3, k3.pos._3, k4.pos._3, H / 6.0);

                output.v._1 = recalculation(output.v._1, k.v._1, k2.v._1, k3.v._1, k4.v._1, H / 6.0);
                output.v._2 = recalculation(output.v._2, k.v._2, k2.v._2, k3.v._2, k4.v._2, H / 6.0);
                output.v._3 = recalculation(output.v._3, k.v._3, k2.v._3, k3.v._3, k4.v._3, H / 6.0);

                R = Math.Sqrt(Math.Pow(output.pos._1, 2) + Math.Pow(output.pos._2, 2) + Math.Pow(output.pos._3, 2));
                t = t + H;// сдвигаемся по времени
                          // пересчет шага H

                double K1 = Math.Pow(k.pos._1 / 6.0 + k2.pos._1 * 2 / 6.0 - 4 * k3.pos._1 / 6.0 + k4.pos._1 / 6.0, 2);
                double K2 = Math.Pow(k.pos._2 / 6.0 + k2.pos._2 * 2 / 6.0 - 4 * k3.pos._2 / 6.0 + k4.pos._2 / 6.0, 2);
                double K3 = Math.Pow(k.pos._3 / 6.0 + k2.pos._3 * 2 / 6.0 - 4 * k3.pos._3 / 6.0 + k4.pos._3 / 6.0, 2);
                double Ecal = Math.Abs(H * Math.Pow(K1 + K2 + K3, 0.5));

                H = H * Math.Pow(EpsPresset / Ecal, 1/3.0);
                qqq++;
                System.Console.WriteLine(H);// отслеживаем шаги
                textFile.WriteLine(H);//вывод шага
                File.WriteLine(R);//вывод гелиоцентрического расстояния
                text.WriteLine(t);// вывод времени

            }
            System.Console.WriteLine(t);//конечный момент времени
            System.Console.WriteLine(qqq);//колличество шагов
            textFile.Close();
            File.Close();
            text.Close();
            return output;
        }

        //персчет координат и скоростей
        static double recalculation(double x, double k1, double k2, double k3, double k4, double h)
        {
            return x + (k1 + 2 * k2 + 2 * k3 + k4) * h;
        }

        // функция которая считает возмущения
        public static double Vozmysheniya(double t, int nom, Vector3 массив_координат_астеройда, DEreader DE405, double jdate)
        {

            int index = 0;//Планета
            bool geleo = true;//гелео=тру  бари=фолс
            jdate = jdate + t;//Дата
            double[] poz = new double[3];//массив куда счит (коор=3) (коор+скор=6)
            double[] Mas = new double[10];//массив где хранятся массы планет
            double del3;
            double r3;
            double G = 2.959122082855911025E-4;
            double EarthDivMoon = 0.813005600000000044E+02;
            Mas[0] = 0.491254745145081187E-10;
            Mas[1] = 0.724345248616270270E-09;
            Mas[2] = 0.899701134671249882E-09 * (1 - 1 / (EarthDivMoon + 1)); // Earth
            Mas[3] = 0.954953510577925806E-10;
            Mas[4] = 0.282534590952422643E-06;
            Mas[5] = 0.845971518568065874E-07;
            Mas[6] = 0.129202491678196939E-07;
            Mas[7] = 0.152435890078427628E-07;
            Mas[8] = 0.218869976542596968E-11;
            Mas[9] = 0.899701134671249882E-09 / (EarthDivMoon + 1);
            double sum = 0;
            while (index < 10)
            {
                DE405.DEreaderGetPlanetPoz(false, jdate, index, geleo, poz);
                r3 = Math.Pow(poz[0] * poz[0] + poz[1] * poz[1] + poz[2] * poz[2], 0.5);
                r3 = r3 * r3 * r3;
                del3 = Math.Pow((массив_координат_астеройда.values[0] - poz[0]) * (массив_координат_астеройда.values[0] - poz[0]) +
                          (массив_координат_астеройда.values[1] - poz[1]) * (массив_координат_астеройда.values[1] - poz[1]) +
                          (массив_координат_астеройда.values[2] - poz[2]) * (массив_координат_астеройда.values[2] - poz[2]), 0.5);
                del3 = del3 * del3 * del3;
                sum += Mas[index] * ((poz[nom] - массив_координат_астеройда.values[nom]) / del3 - poz[nom] / r3);
                index++;
            }
            sum = sum * G;
            return sum;
        }

        static void Main(string[] args)
        {
            DEreader DE405 = new DEreader(405, "C: \\Users\\Amanda\\Desktop\\16002200.405");//бинарный файл
            Coordinates input = new Coordinates();
            double t0, tn; // начальные данные
            double jdate = 2457700.5;//Дата

            //// ввод информации
            System.Console.WriteLine("Введите t0, а затем tn");
            bool success = Double.TryParse(System.Console.ReadLine(), out t0);
            bool success2 = Double.TryParse(System.Console.ReadLine(), out tn);
            if (!success || !success2)
            {
                Console.WriteLine("Не верный ввод данных!");
                return;
            }
            /// считываем координаты из файла
            using (StreamReader fs = new StreamReader((new FileStream(filePath, FileMode.Open))))
            {

                input.pos._1 = double.Parse(fs.ReadLine(), CultureInfo.InvariantCulture);
                input.pos._2 = double.Parse(fs.ReadLine(), CultureInfo.InvariantCulture);
                input.pos._3 = double.Parse(fs.ReadLine(), CultureInfo.InvariantCulture);
                input.v._1 = double.Parse(fs.ReadLine(), CultureInfo.InvariantCulture);
                input.v._2 = double.Parse(fs.ReadLine(), CultureInfo.InvariantCulture);
                input.v._3 = double.Parse(fs.ReadLine(), CultureInfo.InvariantCulture);
                Console.WriteLine(input.pos._1);
                Console.WriteLine(input.pos._2);
                Console.WriteLine(input.pos._3);
                Console.WriteLine(input.v._1);
                Console.WriteLine(input.v._2);
                Console.WriteLine(input.v._3);
            }

            input = calc(input, t0, tn, DE405, jdate); //запускаем программу на счет

            // выводим результат
            using (StreamWriter sw = new StreamWriter(new FileStream(output, FileMode.Create)))
            {
                sw.WriteLine("x=" + input.pos._1);
                sw.WriteLine("y=" + input.pos._2);
                sw.WriteLine("z=" + input.pos._3);
                sw.WriteLine("Vx=" + input.v._1);
                sw.WriteLine("Vy=" + input.v._2);
                sw.WriteLine("Vz=" + input.v._3);
            }

            Console.ReadLine();
        }
    }
}