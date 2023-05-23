using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Xml.Linq;
using Color = System.Drawing.Color;

namespace FDMApp
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            Bitmap bitmap = new Bitmap(Config.ScreenWidth, Config.ScreenHeight);
            double[,] nums = new double[bitmap.Width, bitmap.Height];
            nums.Initialize();
            List<List<Node>> nodes = new List<List<Node>>();
            InitializeComponent();
            CreateGrid(nodes);
            Solve(ref nodes);
            var sl = nodes.SelectMany(q => q);
            double min = sl.Where(t => t.type == 0).Min(v => v.val);
            double max = sl.Where(t => t.type == 0).Max(v => v.val);
            Interp(ref nums, nodes);
            for (int i = 0; i < Config.ScreenWidth; i++)
            {
                for (int j = 0; j < Config.ScreenHeight; j++)
                {
                    if (nums[i, j] < min) nums[i, j] = min;
                    if (nums[i, j] > max) nums[i, j] = max;
                }
            }
            CreateBitmap(ref bitmap, ref nums, min, max);
            FillBlankSpace(ref bitmap);
            PrintToFile(nodes);
            ApplyNodes(ref bitmap, ref nodes);
            this.MainImage.Source = BitmapToImageSource(bitmap);
        }
        static void PrintToFile(List<List<Node>> nodes)
        {
            string path = $"Log{DateTime.Now.Second}.txt";
            string text = "";
            List<List<string>> printNode = new List<List<string>>();
            for (double y = Config.VerticalLength - Config.DeltaY; y > 0; y -= Config.DeltaY)
            {
                for (double x = 0; x < Config.HorizontaLength - Config.DeltaX; x += Config.DeltaX)
                {
                    Node curNode = nodes.GetNodeWith(x, y);
                    if (Math.Sqrt(Math.Pow(curNode.x - x, 2) + Math.Pow(curNode.y - y, 2)) < Config.DeltaX)
                        text += String.Format("{0:f2} ", curNode.val);
                    else text += "     ";

                }
                text += "\n";
            }
            using (FileStream fstream = new FileStream(path, FileMode.OpenOrCreate))
            {
                byte[] input = Encoding.Default.GetBytes(text);
                fstream.Write(input, 0, input.Length);
            }
        }
        static void Solve(ref List<List<Node>> nodes)
        {
            double l = 0;
            double n = 0;
            double k = Config.KVal;
            double dx = Config.DeltaX; double dy = Config.DeltaY;
            int nCountX = (int)(Config.HorizontaLength / dx);
            int nCountY = (int)(Config.VerticalLength / dy);
            double dVal = 0; double aVal = 0; double bVal = 0; double cVal = 0;
            double dt = Config.DeltaTime;
            double x; double y; Node boundNode;
            for (double t = 0; t < Config.Time; t += dt)
            {
                int count = 0;
                List<double> a = new List<double>();
                List<double> b = new List<double>();
                List<double> c = new List<double>();
                List<double> d = new List<double>();
                List<double> ans = new List<double>();
                for (int i = 0; i < nCountY; i++)
                {
                    for (int j = 0; j < nCountX; j++)
                    {
                        dVal = 0; aVal = 0; bVal = 0; cVal = 0;
                        var des = nodes.GetNodeWith(j * dx, i * dy);
                        if (Distance(des, j * dx, i * dy) > dx / 10d) continue;
                        if (des.type > 0) { continue; }
                        x = des.x;
                        y = des.y;
                        count++;
                        boundNode = nodes.GetNodeWith(x - dx, y);
                        n = Math.Abs((boundNode.x - des.x)) / dx;
                        n = n > 0.99d ? 1d : n;
                        n = n < 0.99d ? 0d : n;
                        if (boundNode.type > 0)
                        {
                            boundNode.val = Config.GetBoundaryValue(boundNode, des, dx);
                            if (n > 0 && n < 1)
                                dVal += (dt * k * n / (dx * dx * (n * (n + 1)))) * boundNode.val;
                            else
                                dVal += (dt * k / (dx * dx)) * boundNode.val;
                            aVal = 0;
                        }
                        else
                        {
                            if (n > 0 && n < 1)
                                aVal = -dt * k * n / (dx * dx * (n * (n + 1)));
                            else
                                aVal = -dt * k / (dx * dx);
                        }
                        boundNode = nodes.GetNodeWith(x + dx, y);
                        if (boundNode.type > 0)
                        {
                            boundNode.val = Config.GetBoundaryValue(boundNode, des, dx);
                            if (n > 0 && n < 1)
                                dVal += (dt * k * n / (dx * dx * (n * (n + 1)))) * boundNode.val;
                            else
                                dVal += (dt * k / (dx * dx)) * boundNode.val;
                            cVal = 0;
                        }
                        else
                        {
                            if (n > 0 && n < 1)
                                cVal = -dt * k * n / (dx * dx * (n * (n + 1)));
                            else
                                cVal = -dt * k / (dy * dx);
                        }
                        if (n > 0 && n < 1)
                            bVal = (2 * dt * k * (n + 1) / (dx * dx * (n * (n + 1))) + 1);
                        else
                            bVal = (2 * dt * k / (dx * dx) + 1);
                        a.Add(aVal);
                        b.Add(bVal);
                        c.Add(cVal);
                        d.Add(des.val + dVal);
                    }
                }
                solveMatrix(count, a, b, c, d, ans);
                count = 0;
                for (int i = 0; i < nCountY; i++)
                {
                    for (int j = 0; j < nCountX; j++)
                    {
                        var des = nodes.GetNodeWith(j * dx, i * dy);
                        if (Distance(des, j * dx, i * dy) > dx / 10d) continue;
                        if (des.type > 0) continue;
                        des.val = ans[count];
                        count++;
                    }
                }
                a.Clear(); b.Clear(); c.Clear(); d.Clear(); ans.Clear(); count = 0;
                for (int i = 0; i < nCountX; i++)
                {
                    for (int j = 0; j < nCountY; j++)
                    {
                        dVal = 0; aVal = 0; bVal = 0; cVal = 0;
                        var des = nodes.GetNodeWith(i * dx, j * dy);
                        if (Distance(des, i * dx, j * dy) > dx / 10d) continue;
                        if (des.type > 0) { continue; }
                        x = des.x;
                        y = des.y;
                        count++;
                        boundNode = nodes.GetNodeWith(x, y - dy);
                        l = Math.Abs((boundNode.y - des.y)) / dy;
                        l = l > 0.99d ? 1d : l;
                        l = l < 0.99d ? 0d : l;
                        if (boundNode.type > 0)
                        {
                            boundNode.val = Config.GetBoundaryValue(boundNode, des, dy);
                            if (l > 0 && l < 1)
                                dVal += (dt * k * l / (dy * dy * (l * (l + 1)))) * boundNode.val;
                            else
                                dVal += (dt * k / (dy * dy)) * boundNode.val;
                            aVal = 0;
                        }
                        else
                        {
                            if (l > 0 && l < 1)
                                aVal = -dt * k * l / (dy * dy * (l * (l + 1)));
                            else
                                aVal = -dt * k / (dy * dy);
                        }
                        boundNode = nodes.GetNodeWith(x, y + dy);
                        l = Math.Abs((boundNode.y - des.y)) / dy;
                        l = l > 0.99d ? 1d : l;
                        l = l < 0.99d ? 0d : l;
                        if (boundNode.type > 0)
                        {
                            boundNode.val = Config.GetBoundaryValue(boundNode, des, dy);
                            if (l > 0 && l < 1)
                                dVal += (dt * k * l / (dy * dy * (l * (l + 1)))) * boundNode.val;
                            else
                                dVal += (dt * k / (dy * dy)) * boundNode.val;
                            cVal = 0;
                        }
                        else
                        {
                            cVal = -dt * k / (dy * dy);
                        }
                        if (l > 0 && l < 1)
                            bVal = (2 * dt * k * (l + 1) / (dy * dy * (l * (l + 1))) + 1);
                        else
                            bVal = (2 * dt * k / (dy * dy) + 1);
                        a.Add(aVal);
                        b.Add(bVal);
                        c.Add(cVal);
                        d.Add(des.val + dVal);
                    }
                }
                solveMatrix(count, a, b, c, d, ans);
                count = 0;
                for (int i = 0; i < nCountX; i++)
                {
                    for (int j = 0; j < nCountY; j++)
                    {
                        var des = nodes.GetNodeWith(i * dx, j * dy);
                        if (Distance(des, i * dx, j * dy) > dx / 10d) continue;
                        if (des.type > 0) continue;
                        des.val = ans[count];
                        count++;
                    }
                }
            }
        }
        public static void solveMatrix(int n, List<double> a, List<double> b, List<double> c, List<double> d, List<double> x)
        {
            var y = new List<double>();
            var alpha = new List<double>();
            var betha = new List<double>();
            y.Add(b[0]);
            alpha.Add(-c[0] / y[0]);
            betha.Add(d[0] / y[0]);
            for (int i = 1; i < n - 1; i++)
            {
                y.Add(b[i] + a[i] * alpha[i - 1]);
                alpha.Add(-c[i] / y[i]);
                betha.Add((d[i] - a[i] * betha[i - 1]) / y[i]);
            }
            y.Add(b[n - 1] + a[n - 1] * alpha[n - 2]);
            betha.Add((d[n - 1] - a[n - 1] * betha[n - 2]) / y[n - 1]);
            for (int i = 0; i < n; i++)
            {
                x.Add(0);
            }
            x[n - 1] = betha[n - 1];
            for (int i = n - 2; i >= 0; i--)
            {
                x[i] = (alpha[i] * x[i + 1] + betha[i]);
            }
        }
        static void ApplyNodes(ref Bitmap bitmap, ref List<List<Node>> nodes)
        {
            double dx = Config.HorizontaLength / Config.ScreenWidth;
            double dy = Config.VerticalLength / Config.ScreenHeight;
            foreach (List<Node> nodeList in nodes)
            {
                foreach (Node node in nodeList)
                {
                    if (node.type == BoundType.Inside)
                        bitmap.SetPixel((int)(node.x / dx), (int)(node.y / dy), Color.White);
                    if (node.type == BoundType.IsInBigCircle)
                        bitmap.SetPixel((int)(node.x / dx), (int)(node.y / dy), Color.Black);
                    if (node.type == BoundType.IsInSmallCircle)
                        bitmap.SetPixel((int)(node.x / dx), (int)(node.y / dy), Color.Red);
                    if (node.type == BoundType.InLeftSide)
                        bitmap.SetPixel((int)(node.x / dx), (int)(node.y / dy), Color.Yellow);
                    if (node.type == BoundType.InRightSide)
                        bitmap.SetPixel((int)(node.x / dx), (int)(node.y / dy), Color.Yellow);
                    if (node.type == BoundType.InUpperSide)
                        bitmap.SetPixel((int)(node.x / dx), (int)(node.y / dy), Color.Blue);
                    if (node.type == BoundType.InLowerSide)
                        bitmap.SetPixel((int)(node.x / dx), (int)(node.y / dy), Color.Blue);
                }
            }
        }
        static void FillBlankSpace(ref Bitmap bitmap)
        {
            double dx = Config.HorizontaLength / Config.ScreenWidth;
            double dy = Config.VerticalLength / Config.ScreenHeight;
            for (int x = 0; x < bitmap.Width; ++x)
            {
                for (int y = 0; y < bitmap.Height; ++y)
                {
                    if (NodePos(x * dx, y * dy) < 0) bitmap.SetPixel(x, y, Color.White);
                }
            }
        }
        static void CreateGrid(List<List<Node>> nodes)
        {
            double dx = Config.DeltaX; double dy = Config.DeltaY;
            double maxX = Config.HorizontaLength; double maxY = Config.VerticalLength;
            double x = 0;
            BoundType boundType;
            for (int i = 0; x < maxX; i++, x += dx)
            {
                double y = 0;
                nodes.Add(new List<Node>());
                for (int j = 0; y < maxY; j++, y += dy)
                {
                    if ((boundType = NodePos(x, y)) >= 0)
                    {
                        if (boundType == BoundType.Inside)
                            nodes[i].Add(new Node(x, y, Config.InitPlateT, boundType));
                        if (boundType == BoundType.IsInBigCircle)
                            nodes[i].Add(new Node(x, y, Config.InitBigCircleT, boundType));
                        if (boundType == BoundType.IsInSmallCircle)
                            nodes[i].Add(new Node(x, y, Config.InitSmallCircleT, boundType));
                        if (boundType == BoundType.InLeftSide)
                            nodes[i].Add(new Node(x, y, Config.InitLeftT, boundType));
                        if (boundType == BoundType.InRightSide)
                            nodes[i].Add(new Node(x, y, Config.InitRightSide, boundType));
                        if (boundType == BoundType.InUpperSide)
                            nodes[i].Add(new Node(x, y, Config.InitUpperSide, boundType));
                        if (boundType == BoundType.InLowerSide)
                            nodes[i].Add(new Node(x, y, Config.InitLowerSide, boundType));
                    }
                }
            }
            nodes.GetNodeWith(0d, 6d).val = Config.InitLeftT;
        }
        static BoundType NodePos(double x, double y)
        {
            double r1 = Config.SmallCircleRad;
            double r2 = Config.BigCircleRad;
            double dx = Config.DeltaX;
            double dy = Config.DeltaY;
            if (Math.Pow(x - 2, 2) + Math.Pow(y - 4, 2) < Math.Pow(r1, 2))
            {
                return BoundType.None;
            }
            if (Math.Pow((x + dx) - 2, 2) + Math.Pow((y) - 4, 2) < Math.Pow(r1, 2) || Math.Pow((x) - 2, 2) + Math.Pow((y + dy) - 4, 2) < Math.Pow(r1, 2)) return BoundType.IsInSmallCircle;
            if (Math.Pow((x - dx) - 2, 2) + Math.Pow((y) - 4, 2) < Math.Pow(r1, 2) || Math.Pow((x) - 2, 2) + Math.Pow((y - dy) - 4, 2) < Math.Pow(r1, 2)) return BoundType.IsInSmallCircle;
            if (Math.Pow(x - (Config.HorizontaLength - r2), 2) + Math.Pow(y - (Config.VerticalLength - r2), 2) < Math.Pow(r2, 2))
            {
                if (!(y < (Config.VerticalLength - r2) && x < Config.HorizontaLength) && !(x < (Config.HorizontaLength - r2) && y < Config.VerticalLength))
                {
                    if (Math.Pow((x) - (Config.HorizontaLength - r2), 2) + Math.Pow((y + dy) - (Config.VerticalLength - r2), 2) > Math.Pow(r2, 2)) return BoundType.IsInBigCircle;
                    if (Math.Pow((x + dx) - (Config.HorizontaLength - r2), 2) + Math.Pow((y) - (Config.VerticalLength - r2), 2) > Math.Pow(r2, 2)) return BoundType.IsInBigCircle;
                    return BoundType.Inside;
                }
            }
            if (x - dx < 0 && y + dy > Config.VerticalLength) return BoundType.InUpperSide;
            if (x - dx < 0)
                return BoundType.InLeftSide;
            if (y - dy < 0) return BoundType.InLowerSide;
            if (y < (Config.VerticalLength - r2) && x < Config.HorizontaLength)
            {
                if (x + dx >= Config.HorizontaLength) return BoundType.InRightSide;
                return BoundType.Inside;
            }
            if (x < (Config.HorizontaLength - r2) && y < Config.VerticalLength)
            {
                if (y + dy > Config.VerticalLength) return BoundType.InUpperSide;
                return BoundType.Inside;
            }
            return BoundType.None;
        }
        static double Distance(Node n, double x, double y)
        {
            return (n.x - x) * (n.x - x) + (n.y - y) * (n.y - y);
        }
        static void Interp(ref double[,] nums, List<List<Node>> nodes)
        {
            List<Node> nodes1 = new List<Node>();
            double x = 0; double y = 0;
            double dx = Config.HorizontaLength / Config.ScreenWidth;
            double dy = Config.VerticalLength / Config.ScreenHeight;
            foreach (List<Node> i in nodes)
            {
                foreach (Node j in i)
                {
                    nodes1.Add(j);
                }
            }
            for (int i = 0; i < nums.GetLength(0); i++, x += dx)
            {
                y = 0;
                for (int j = 0; j < nums.GetLength(1); j++, y += dy)
                {
                    var sorted = nodes1
                        .OrderBy(n => Distance(n, x, y)).Take(3).ToList();
                    double d = (sorted[1].y - sorted[2].y) * (sorted[0].x - sorted[2].x) + (sorted[2].x - sorted[1].x) * (sorted[0].y - sorted[2].y);
                    if (NodePos(x, y) < 0 || d == 0) { nums[i, j] = sorted[0].val; continue; }
                    double t1 = ((sorted[1].y - sorted[2].y) * (x - sorted[2].x) + (sorted[2].x - sorted[1].x) * (y - sorted[2].y)) / d;
                    double t2 = ((sorted[2].y - sorted[0].y) * (x - sorted[2].x) + (sorted[0].x - sorted[2].x) * (y - sorted[2].y)) / d;
                    double t3 = 1 - t1 - t2;
                    nums[i, j] = t1 * sorted[0].val + t2 * sorted[1].val + t3 * sorted[2].val;
                }
            }
        }

        static void CreateBitmap(ref Bitmap bitmap, ref double[,] nums, double minT, double maxT)
        {
            for (int x = 0; x < bitmap.Width; ++x)
            {
                for (int y = 0; y < bitmap.Height; ++y)
                {
                    double h = (240d / 360d) - (nums[x, y] - minT) / (maxT - minT) * 240d / 360d;
                    var rgb = ColorUtils.HSVToRGB(h, 1, 1);
                    bitmap.SetPixel(x, y, rgb);
                }
            }
        }
        private BitmapImage BitmapToImageSource(System.Drawing.Bitmap bitmap)
        {
            using (MemoryStream memory = new MemoryStream())
            {
                bitmap.Save(memory, System.Drawing.Imaging.ImageFormat.Bmp);
                memory.Position = 0;
                BitmapImage bitmapimage = new BitmapImage();
                bitmapimage.BeginInit();
                bitmapimage.StreamSource = memory;
                bitmapimage.CacheOption = BitmapCacheOption.OnLoad;
                bitmapimage.EndInit();
                return bitmapimage;
            }
        }
    }
}
