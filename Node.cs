using System.Collections.Generic;
using System;
using System.Linq;

namespace FDMApp
{
    class Node
    {
        public double x { get; set; }
        public double y { get; set; }
        public double val { get; set; }
        public bool IsBound { get; set; } = false;
        public BoundType type { get; set; }
        public Node(double x, double y, double val, BoundType bt)
        {
            this.x = x;
            this.y = y;
            this.val = val;
            this.type = bt;
        }
    }
    public enum BoundType : int
    {
        None = -1,
        Inside = 0,
        IsInBigCircle = 1,
        IsInSmallCircle = 2,
        InRightSide = 3,
        InLeftSide = 4,
        InUpperSide = 5,
        InLowerSide = 6
    }
    static class ExtensionNode
    {
        public static Node GetNodeWith(this List<List<Node>> nodes, double x, double y)
        {
            double dist = int.MaxValue;
            Node val = null;
            for (int i = 0; i < nodes.Count(); i++)
            {
                for (int j = 0; j < nodes[i].Count(); j++)
                {
                    //if (nodes[i][j].type == 0) continue;
                    double d = Distance(nodes[i][j], x, y);
                    if (d < dist)
                    {
                        val = nodes[i][j];
                        dist = d;
                    }
                }
            }
            return val;
        }
        static double Distance(Node n, double x, double y)
        {
            return (n.x - x) * (n.x - x) + (n.y - y) * (n.y - y);
        }
    }
}

