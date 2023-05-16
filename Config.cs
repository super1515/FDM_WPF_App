using System;

namespace FDMApp
{
    static internal class Config
    {
        public const double InitPlateT = 10;
        public const double InitUpperSide = 10;
        public const double InitRightSide = 10;
        public const double InitLowerSide = 10;
        public const double InitLeftT = 20;
        public const double InitBigCircleT = 100;
        public const double InitSmallCircleT = 10;

        public const double BigCircleRad = 3d;
        public const double SmallCircleRad = 1d;

        public const double HorizontaLength = 8;
        public const double VerticalLength = 6;

        public const double KVal = 50d;

        public const double DeltaX = 0.2d;
        public const double DeltaY = 0.2d;
        public const double DeltaTime = 0.1d;
        public const double Time = 25d;

        public const int ScreenWidth = 400;
        public const int ScreenHeight = 300;
        public static double GetBoundaryValue(Node boundNode, Node curNode, double delta)
        {
            switch (boundNode.type)
            {
                case BoundType.IsInBigCircle:
                    return InitBigCircleT;
                case BoundType.IsInSmallCircle:
                    return InitSmallCircleT;
                case BoundType.InLeftSide:
                    return InitLeftT;
                case BoundType.InRightSide:
                    return -20 * delta + curNode.val;
                case BoundType.InUpperSide:
                    return -delta * curNode.val + curNode.val;
                case BoundType.InLowerSide:
                    return -delta * curNode.val + curNode.val;
            }
            throw new ArgumentException();
        }
    }
}
