using System;
using System.Drawing;

namespace FDMApp
{
    internal class ColorUtils
    {

        public static Color HSVToRGB(double H, double S, double V)
        {
            double R = 0, G = 0, B = 0;
            if (H == 1.0)
            {
                H = 0.0;
            }

            double step = 1.0 / 6.0;
            double vh = H / step;

            int i = (int)System.Math.Floor(vh);

            double f = vh - i;
            double p = V * (1.0 - S);
            double q = V * (1.0 - (S * f));
            double t = V * (1.0 - (S * (1.0 - f)));

            switch (i)
            {
                case 0:
                    {
                        R = V;
                        G = t;
                        B = p;
                        break;
                    }
                case 1:
                    {
                        R = q;
                        G = V;
                        B = p;
                        break;
                    }
                case 2:
                    {
                        R = p;
                        G = V;
                        B = t;
                        break;
                    }
                case 3:
                    {
                        R = p;
                        G = q;
                        B = V;
                        break;
                    }
                case 4:
                    {
                        R = t;
                        G = p;
                        B = V;
                        break;
                    }
                case 5:
                    {
                        R = V;
                        G = p;
                        B = q;
                        break;
                    }
                default:
                    {
                        // not possible - if we get here it is an internal error
                        throw new ArgumentException();
                    }
            }
            return Color.FromArgb((int)(R * 255), (int)(G * 255), (int)(B * 255));
        }
    }
}

