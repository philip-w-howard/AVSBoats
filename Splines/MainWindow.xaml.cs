using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Diagnostics;

namespace Splines
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private Point[] pointList = new Point[] { };
        private int numDots = 0;

        private void Canvas_MouseDown(object sender, MouseButtonEventArgs e)
        {
            int size = 10;
            Point location = e.GetPosition(drawingCanvas);
            Debug.WriteLine(location);
            pointList[numDots] = location;
            numDots++;

            Ellipse _circle = new Ellipse
            {
                Width = size,
                Height = size,
                Stroke = Brushes.Black,
                StrokeThickness = 1
            };
            Canvas.SetLeft(_circle, e.GetPosition(drawingCanvas).X - size/2);
            Canvas.SetTop(_circle, e.GetPosition(drawingCanvas).Y - size/2);
            drawingCanvas.Children.Add(_circle);
        }

        private void Draw_Click(object sender, RoutedEventArgs e)
        {
            //Pen pen = new Pen(Brushes.Black);
            //drawingCanvas.Drawing.DrawCurve(pen, pointList, 0.5f);
        }
    }
}
