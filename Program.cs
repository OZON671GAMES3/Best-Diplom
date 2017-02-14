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


//-----------------------------------------------------------------------------------------------Старое
    class Program
    {
       public double[] Kut(double[] components, double tn, double eps,double h0, DEreader DE405,
           double jdate)// метод рунгекутта
        {
            double t = 0;
            double h = h0;
            double[] fake = new double[3]; fake[0] = 0; fake[1] = 0; fake[2] = 0;
            while (t < tn)
            {
                if (t + h > tn)
                { h = tn - t; };
                double[] k1 = fcn(components, h/2, t, fake,DE405,jdate);
                double[] k2 = fcn(components, h/2, t,k1, DE405, jdate);
                double[] k3 = fcn(components, h/2, t,k2, DE405, jdate);
                double[] k4 = fcn(components, h, t,k3, DE405, jdate);
                components = newpos(components, k1, k2, k3, k4, h / 6);
                t += h;
                // пересчет шага H

                double K1 = Math.Pow(k1[1] / 6 + k2[1] * 2 / 6 - 4 * k3[1] / 6 + k4[1] / 6, 2);
                double K2 = Math.Pow(k1[2] / 6 + k2[2] * 2 / 6 - 4 * k3[2] / 6 + k4[2] / 6, 2);
                double K3 = Math.Pow(k1[3] / 6 + k2[3] * 2 / 6 - 4 * k3[3] / 6 + k4[3] / 6, 2);
                double Ecal = Math.Abs(h * Math.Pow(K1 + K2 + K3, 0.5));

                h = h * Math.Pow(eps / Ecal, 0.33333);
            }
            return components;
        }

        //функция правых частей
        public double[] fcn(double[] components, double h, double t, double[] fake,DEreader DE405,double jdate)
        {
            double G = 2.959122082855911025E-4;
            //подсчет радиуса
            double R = Math.Sqrt(Math.Pow((components[0] + fake[0] * h), 2) + Math.Pow((components[1] + fake[1] * h), 2) +
                Math.Pow((components[2] + fake[2] * h), 2));
            double[] vektor = new double[6];
            //часть связаная с возмущениями
            bool geleo = true;//гелео=тру  бари=фолс
            jdate = jdate + t;
            double[] poz = new double[3];//массив куда счит (коор=3) (коор+скор=6)
            double[] Mas = new double[10];//массив где хранятся массы планет
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
            //пересчет компонент
            for (int i = 0; i < 3; i++)
            {
            vektor[i] = components[i + 3];// определяем позицию
            double sum = 0;
            int index = 0;
            double r3;// расстояние от центра в кубе xj/r^3
            double del3;// треугольник в кубе, разница позиций планеты и астеройда
            while (index < 10) //считаем возмущение
            {
               DE405.DEreaderGetPlanetPoz(false, jdate, index, geleo, poz);// cчитываем координаты планеты на момент времени
               r3 = Math.Pow(poz[0] * poz[0] + poz[1] * poz[1] + poz[2] * poz[2], 0.5);
               r3 = r3 * r3 * r3;
               del3 = Math.Pow((components[0] - poz[0]) * (components[0] - poz[0]) +
                     (components[1] - poz[1]) * (components[1] - poz[1]) +
                     (components[2] - poz[2]) * (components[2] - poz[2]), 0.5);
               del3 = del3 * del3 * del3;
               sum += Mas[index] * ((poz[i] - components[i]) / del3 - poz[i] / r3);
               index++;
            }
                vektor[i+3]= (-G * components[i]) / (R * R * R)+sum;// скорость +возмущение
            }
            return vektor;
        }
        //персчет координат и скоростей
        public double[] newpos(double[] components, double[] k1, double[] k2, double[] k3, double[] k4, double h)
        {
            double[] vektor = new double[6];
            for(int i=0; i<6; i++)
            {
                vektor[i] = components[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * h;
            }
            return vektor;
        }


        static void Main(string[] args)
        {
            DEreader DE405 = new DEreader(405, "C: \\Users\\Amanda\\Desktop\\16002200.405");//бинарный файл
            double jdate = 2457700.5;//Дата
            double[] components = new double[6];
            double[] otvet= new double[6];
            double t,h0,tn,eps;
            h0 = 0.001;
            eps = 1E-12;
            tn = 1680;
            components[0]= 2.43217016648123;
            components[1] = 1.50652966672052;
            components[2] = 0.214850898091387;
            components[3] = -0.00551948749935559;
            components[4] = 0.00701865011326745;
            components[5] = 0.00443325211708482;




            otvet =Kut(components, tn, eps, h0, DE405, jdate);

            Console.ReadLine();
        }
    }
}