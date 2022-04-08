using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
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

namespace Biometria_1
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        //public readonly static string ImageFilePath = "../../../A.jpg";
        public readonly static string ImageFilePath = "../../../testowy.jpg";
        //public readonly static string ImageFilePath = "../../../Tree.jpg";
        public readonly static string ResultFilePath = "../../../Result.jpg";

        public Bitmap BitmapSrc;
        public Bitmap BitmapResult;
        public MainWindow()
        {
            BitmapSrc = new Bitmap(ImageFilePath);
            InitializeComponent();

            Threshold.Minimum = byte.MinValue;
            Threshold.Maximum = byte.MaxValue;

            StretchingMin.Minimum = 0;
            StretchingMin.Maximum = 254;
            StretchingMax.Minimum = 1;
            StretchingMax.Maximum = 255;

            UpdateImage();
        }

        public enum HistogramMode
        {
            Red,
            Green,
            Blue,
            Average
        }

        public Bitmap Histogram(Bitmap bitmap, HistogramMode mode)
        {
            var data = bitmap.LockBits(
                new System.Drawing.Rectangle(new System.Drawing.Point(0, 0), bitmap.Size),
                System.Drawing.Imaging.ImageLockMode.ReadWrite,
                System.Drawing.Imaging.PixelFormat.Format24bppRgb);

            var bitmapData = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapData, 0, bitmapData.Length);

            int[] histogramValues = GetHistogram(mode, bitmapData);

            double maxValue = histogramValues.Max();

            for (int i = 0; i < histogramValues.Length; i++)
            {
                histogramValues[i] = (int)(histogramValues[i] / maxValue * data.Height);
            }

            bitmapData = new byte[bitmapData.Length];
            for (int i = 0; i < bitmapData.Length; i++)
            {
                bitmapData[i] = byte.MaxValue;
            }


            for (int i = 0; i < histogramValues.Length; i++)
            {
                for (int j = 0; j < histogramValues[i]; j++)
                {
                    int index = i * 3 + (data.Height - 1 - j) * data.Stride + data.Width / byte.MaxValue * i;

                    bitmapData[index] = bitmapData[index + 1] = bitmapData[index + 2] = 0;
                }
            }

            Marshal.Copy(bitmapData, 0, data.Scan0, bitmapData.Length);

            bitmap.UnlockBits(data);

            return bitmap;
        }

        private int[] GetHistogram(HistogramMode mode, byte[] bitmapData)
        {
            int[] histogramValues = new int[byte.MaxValue + 1];
            if (mode == HistogramMode.Average)
            {
                for (int i = 0; i < bitmapData.Length; i += 3)
                {
                    int value = (bitmapData[i] + bitmapData[i + 1] + bitmapData[i + 2]) / 3;
                    histogramValues[value]++;
                }
            }
            else if (mode == HistogramMode.Blue)
            {
                for (int i = 0; i < bitmapData.Length; i += 3)
                {
                    histogramValues[bitmapData[i]]++;
                }
            }
            else if (mode == HistogramMode.Green)
            {
                for (int i = 1; i < bitmapData.Length; i += 3)
                {
                    histogramValues[bitmapData[i]]++;
                }
            }
            else if (mode == HistogramMode.Red)
            {
                for (int i = 2; i < bitmapData.Length; i += 3)
                {
                    histogramValues[bitmapData[i]]++;
                }
            }
            else
            {
                throw new Exception("Wrong histogram mode");
            }
            return histogramValues;
        }


        
        public Bitmap BinaryThreshold(Bitmap bitmap)
        {
            return BinaryThreshold(bitmap, (byte)Threshold.Value);
        }

        public Bitmap BinaryThreshold(Bitmap bitmap, byte threshold, bool average)
        {
            return BinaryThreshold(bitmap, threshold, false, false, false, average);
        }

        public Bitmap BinaryThreshold(Bitmap bitmap, byte threshold)
        {
            return BinaryThreshold(bitmap,
                threshold,
                RedCheckBox.IsChecked.Value,
                GreenCheckBox.IsChecked.Value,
                BlueCheckBox.IsChecked.Value,
                AverageCheckBox.IsChecked.Value);
        }

        public Bitmap BinaryThreshold(Bitmap bitmap, byte threshold, bool red, bool green, bool blue, bool average)
        {
            var data = bitmap.LockBits(
                new System.Drawing.Rectangle(new System.Drawing.Point(0, 0), bitmap.Size),
                System.Drawing.Imaging.ImageLockMode.ReadWrite,
                System.Drawing.Imaging.PixelFormat.Format24bppRgb);

            var bitmapData = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapData, 0, bitmapData.Length);

            for (int i = 0; i < bitmapData.Length; i += 3)
            {
                byte b = blue || average ? bitmapData[i] : byte.MinValue;
                byte g = green || average ? bitmapData[i + 1] : byte.MinValue;
                byte r = red || average ? bitmapData[i + 2] : byte.MinValue;

                byte result = (byte)((r + g + b) / 3);
                if (average)
                {
                    bitmapData[i] = result > threshold ? byte.MaxValue : byte.MinValue;
                    bitmapData[i + 1] = result > threshold ? byte.MaxValue : byte.MinValue;
                    bitmapData[i + 2] = result > threshold ? byte.MaxValue : byte.MinValue;
                }
                else
                {
                    bitmapData[i] = b > threshold ? byte.MaxValue : byte.MinValue;
                    bitmapData[i + 1] = g > threshold ? byte.MaxValue : byte.MinValue;
                    bitmapData[i + 2] = r > threshold ? byte.MaxValue : byte.MinValue;
                }
            }

            Marshal.Copy(bitmapData, 0, data.Scan0, bitmapData.Length);

            bitmap.UnlockBits(data);

            return bitmap;
        }

        public BitmapImage ConvertToImage(Bitmap src)
        {
            MemoryStream ms = new MemoryStream();
            ((System.Drawing.Bitmap)src).Save(ms, System.Drawing.Imaging.ImageFormat.Bmp);
            BitmapImage image = new BitmapImage();
            image.BeginInit();
            ms.Seek(0, SeekOrigin.Begin);
            image.StreamSource = ms;
            image.EndInit();
            return image;
        }

        private void CheckHistogram()
        {
            if (HistogramCheckBox.IsChecked.Value)
            {
                HistogramMode histogram;
                ResultMessage.Content = "Histogram: ";
                if (AverageCheckBox.IsChecked.Value)
                {
                    ResultMessage.Content += "Average";
                    histogram = HistogramMode.Average;
                }
                else if (RedCheckBox.IsChecked.Value)
                {
                    ResultMessage.Content += "Red";
                    histogram = HistogramMode.Red;
                }
                else if (GreenCheckBox.IsChecked.Value)
                {
                    ResultMessage.Content += "Green";
                    histogram = HistogramMode.Green;
                }
                else if (BlueCheckBox.IsChecked.Value)
                {
                    ResultMessage.Content += "Blue";
                    histogram = HistogramMode.Blue;
                }
                else
                {
                    ResultMessage.Content += "Average";
                    histogram = HistogramMode.Average;
                }
                BitmapResult = new Bitmap(BitmapSrc);
                if (StretchingCheckBox.IsChecked.Value)
                    BitmapResult = HistogramStretching(BitmapResult);
                if (EqCheckBox.IsChecked.Value)
                    BitmapResult = HistogramEqual(BitmapResult);
                ImageResult.Source = ConvertToImage(Histogram(BitmapResult, histogram));
            }
            else
            {
                UpdateImage();
            }
        }

        private void UpdateImage()
        {
            BitmapResult = new Bitmap(BitmapSrc);
            if (StretchingCheckBox.IsChecked.Value)
                BitmapResult = HistogramStretching(BitmapResult);

            if (EqCheckBox.IsChecked.Value)
                BitmapResult = HistogramEqual(BitmapResult);

            if (OtsuCheckBox.IsChecked.Value)
                BitmapResult = OtsuMethod(BitmapResult);

            if (NiblackCheckBox.IsChecked.Value)
                BitmapResult = NiBlack(BitmapResult, 7);

            if (SauronCheckBox.IsChecked.Value)
                BitmapResult = Sauvola(BitmapResult, 7);

            if (SkyWalkerCheckBox.IsChecked.Value)
                BitmapResult = Phansalkar(BitmapResult, 7);

            if (BernsenCheckBox.IsChecked.Value)
                BitmapResult = Bernsen(BitmapResult, 7);

            if (RedCheckBox.IsChecked.Value || GreenCheckBox.IsChecked.Value || BlueCheckBox.IsChecked.Value || AverageCheckBox.IsChecked.Value)
            {
                ImageResult.Source = ConvertToImage(BinaryThreshold(BitmapResult));
                ResultMessage.Content = "Binary image: ";
                if (AverageCheckBox.IsChecked.Value)
                    ResultMessage.Content += "Average ";
                else
                {
                    if (RedCheckBox.IsChecked.Value)
                        ResultMessage.Content += "Red ";
                    if (GreenCheckBox.IsChecked.Value)
                        ResultMessage.Content += "Green ";
                    if (BlueCheckBox.IsChecked.Value)
                        ResultMessage.Content += "Blue ";
                }
                ResultMessage.Content += "Threshold: " + (byte)Threshold.Value;
            }
            else
            {
                ImageResult.Source = ConvertToImage(BitmapResult);
                ResultMessage.Content = "Normal image";
            }
        }

        private byte[] LockBitmap(Bitmap bitmap, ref System.Drawing.Imaging.BitmapData data)
        {
            data = bitmap.LockBits(
                new System.Drawing.Rectangle(new System.Drawing.Point(0, 0), bitmap.Size),
                System.Drawing.Imaging.ImageLockMode.ReadWrite,
                System.Drawing.Imaging.PixelFormat.Format24bppRgb);

            var bitmapData = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapData, 0, bitmapData.Length);

            return bitmapData;
        }

        
        private void UpdateImage(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            CheckHistogram();
        }

        private void UpdateImage(object sender, RoutedEventArgs e)
        {
            CheckHistogram();
        }

        private void SaveImage(object sender, RoutedEventArgs e)
        {
            BitmapResult.Save(ResultFilePath);
        }

        private void CheckHistogram(object sender, RoutedEventArgs e)
        {
            CheckHistogram();
        }

        private int[] calculateLUTequal(int[] values, int size)
        {
            double minValue = values.Min();
            int[] output = new int[256];
            double Dn = 0;
            for (int i = 0; i < 256; i++)
            {
                Dn += values[i];
                output[i] = (int)(((Dn - minValue) / (size - minValue)) * 255.0);
            }

            return output;
        }

        private Bitmap HistogramEqual(Bitmap bitmap)
        {
            int[] red = new int[256];
            int[] green = new int[256];
            int[] blue = new int[256];

            for (int x = 0; x < bitmap.Width; x++)
            {
                for (int y = 0; y < bitmap.Height; y++)
                {
                    System.Drawing.Color color = bitmap.GetPixel(x, y);
                    red[color.R]++;
                    green[color.G]++;
                    blue[color.B]++;
                }
            }

            int[] LUTred = calculateLUTequal(red, bitmap.Width * bitmap.Height);
            int[] LUTgreen = calculateLUTequal(green, bitmap.Width * bitmap.Height);
            int[] LUTblue = calculateLUTequal(blue, bitmap.Width * bitmap.Height);
            Bitmap newBitmap = new Bitmap(bitmap.Width, bitmap.Height, bitmap.PixelFormat);


            for (int x = 0; x < bitmap.Width; x++)
            {
                for (int y = 0; y < bitmap.Height; y++)
                {
                    System.Drawing.Color pixel = bitmap.GetPixel(x, y);
                    System.Drawing.Color newpixel = System.Drawing.Color.FromArgb(LUTred[pixel.R], LUTgreen[pixel.G], LUTblue[pixel.B]);
                    newBitmap.SetPixel(x, y, newpixel);
                }
            }
            return newBitmap;
        }
        public Bitmap HistogramEqualization(Bitmap bitmap)
        {
            var data = bitmap.LockBits(
                new System.Drawing.Rectangle(new System.Drawing.Point(0, 0), bitmap.Size),
                System.Drawing.Imaging.ImageLockMode.ReadWrite,
                System.Drawing.Imaging.PixelFormat.Format24bppRgb);

            var bitmapData = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapData, 0, bitmapData.Length);
            bitmap.UnlockBits(data);

            double sum = 0;
            var minValue = bitmapData.Min();

            int[] LUT = new int[bitmapData.Length];

            for (int i = 0; i < bitmapData.Length; i++)
            {
                sum += bitmapData[i];
                LUT[i] = (byte)((float)(sum - minValue) / (bitmapData.Length - minValue) * 255);
            }

            Bitmap newBitmap = new Bitmap(bitmap.Width, bitmap.Height, bitmap.PixelFormat);
            for (int x = 0; x < bitmap.Width; x++)
            {
                for (int y = 0; y < bitmap.Height; y++)
                {
                    System.Drawing.Color pixel = bitmap.GetPixel(x, y);
                    System.Drawing.Color newpixel = System.Drawing.Color.FromArgb(LUT[pixel.B], LUT[pixel.B], LUT[pixel.B]);
                    newBitmap.SetPixel(x, y, newpixel);
                }
            }

            //Marshal.Copy(bitmapData, 0, data.Scan0, bitmapData.Length);

            //bitmap.UnlockBits(data);

            return newBitmap;
        }

        private Bitmap OtsuMethod(Bitmap bitmap)
        {
            System.Drawing.Imaging.BitmapData data = null;
            byte[] bitmapData = LockBitmap(bitmap, ref data);
            int[] histogramValues = GetHistogram(HistogramMode.Average, bitmapData);

            float weightedSumMax = 0;
            for (int i = 0; i < 256; i++)
            {
                weightedSumMax += i * histogramValues[i];
            }

            int total = data.Height * data.Width;
            int sumBefore = 0;
            int threshold = 0;
            float weightedSumBefore = 0;
            float maxvalue = 0;

            for (int i = 0; i < 256; i++)
            {
                sumBefore += histogramValues[i];
                if (sumBefore <= 0) continue;

                int sumAfter = total - sumBefore;
                if (sumAfter <= 0) break;

                weightedSumBefore += (float)(i * histogramValues[i]);

                float mB = weightedSumBefore / sumBefore;
                float mF = (weightedSumMax - weightedSumBefore) / sumAfter;

                float varBetween = (float)sumBefore * sumAfter * (mB - mF) * (mB - mF);
                if (varBetween > maxvalue)
                {
                    maxvalue = varBetween;
                    threshold = i;
                }
            }
            bitmap.UnlockBits(data);
            return BinaryThreshold(bitmap, (byte)threshold, true);
        }

        public Bitmap HistogramStretching(Bitmap bitmap)
        {
            return HistogramStretching(bitmap, (int)StretchingMin.Value, (int)StretchingMax.Value);
        }

        public Bitmap HistogramStretching(Bitmap bitmap, int minValue, int maxValue)
        {
            System.Drawing.Imaging.BitmapData data = null;
            byte[] bitmapData = LockBitmap(bitmap, ref data);

            for (int i = 0; i < bitmapData.Length; i++)
            {
                int diff = (bitmapData[i] - minValue >= 0) ? (bitmapData[i] - minValue) : 0;
                bitmapData[i] = (byte)((float)diff / (maxValue - minValue) * 255);
            }

            Marshal.Copy(bitmapData, 0, data.Scan0, bitmapData.Length);

            bitmap.UnlockBits(data);

            return bitmap;
        }

        public Bitmap NiBlack(Bitmap bitmap, int w = 2, float k = -0.2f)
        {
            System.Drawing.Imaging.BitmapData data = null;
            byte[] bitmapDataIn = LockBitmap(bitmap, ref data);
            byte[] bitmapDataout = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapDataout, 0, bitmapDataout.Length);

            int dy = data.Height, dx = data.Stride;
            //int imgSize = data.Height * data.Width;
            //imgN = copy(img);

            // Calculate the radius of the neighbourhood
            //int w = (n - 1) / 2;

            // Process the image
            for (int i = w + 1; i < dx - w; i++)
            {
                for (int j = w + 1; j < dy - w; j++)
                {
                    List<double> neighbours = new List<double>();
                    // Extract the neighbourhood area
                    for (int x = i - w; x < i + w; x++)
                    {
                        for (int y = j - w; y < j + w; y++)
                        {
                            float bbb = bitmapDataIn[x + y * data.Stride] + bitmapDataIn[x + y * data.Stride + 1] + bitmapDataIn[x + y * data.Stride + 2];
                            bbb /= 3;
                            neighbours.Add(bbb);
                        }
                    }
                    //block = bitmapData[i - w:i + w, j - w:j + w];

                    // Calculate the mean and standard deviation of the neighbourhood region
                    float wBmn = (float)Median(neighbours);
                    float wBstd = (float)standardDeviation(neighbours);

                    // Calculate the threshold value
                    float wBTH = (wBmn + k * wBstd);

                    // Threshold the pixel
                    float aaa = bitmapDataIn[i + j * data.Stride] + 
                        bitmapDataIn[i + j * data.Stride + 1] + 
                        bitmapDataIn[i + j * data.Stride + 2];
                    aaa /= 3;
                    bitmapDataout[i + j * data.Stride] =
                        bitmapDataout[i + j * data.Stride + 1] =
                        bitmapDataout[i + j * data.Stride + 2] = 
                        aaa < wBTH ? byte.MinValue : byte.MaxValue;
                }
            }

            Marshal.Copy(bitmapDataout, 0, data.Scan0, bitmapDataout.Length);
            bitmap.UnlockBits(data);

            return bitmap;
        }

        public double Median(List<double> numbers)
        {
            if (numbers.Count == 0)
                return 0;

            numbers = numbers.OrderBy(n => n).ToList();

            var halfIndex = numbers.Count() / 2;

            if (numbers.Count() % 2 == 0)
                return (numbers[halfIndex] + numbers[halfIndex - 1]) / 2.0;

            return numbers[halfIndex];
        }

        static double standardDeviation(IEnumerable<double> sequence)
        {
            double result = 0;

            if (sequence.Any())
            {
                double average = sequence.Average();
                double sum = sequence.Sum(d => Math.Pow(d - average, 2));
                result = Math.Sqrt((sum) / sequence.Count());
            }
            return result;
        }


        public Bitmap Sauvola(Bitmap bitmap, int w = 2 ,float k = 0.2f, int R = 125)
        {
            System.Drawing.Imaging.BitmapData data = null;
            byte[] bitmapDataIn = LockBitmap(bitmap, ref data);
            byte[] bitmapDataout = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapDataout, 0, bitmapDataout.Length);

            int dy = data.Height, dx = data.Stride;
            //int imgSize = data.Height * data.Width;
            //imgN = copy(img);

            // Calculate the radius of the neighbourhood
            //int w = (n - 1) / 2;

            // Process the image
            for (int i = w + 1; i < dx - w; i++)
            {
                for (int j = w + 1; j < dy - w; j++)
                {
                    List<double> neighbours = new List<double>();
                    // Extract the neighbourhood area
                    for (int x = i - w; x < i + w; x++)
                    {
                        for (int y = j - w; y < j + w; y++)
                        {
                            float bbb = bitmapDataIn[x + y * data.Stride] + bitmapDataIn[x + y * data.Stride + 1] + bitmapDataIn[x + y * data.Stride + 2];
                            bbb /= 3;
                            neighbours.Add(bbb);
                        }
                    }
                    //block = bitmapData[i - w:i + w, j - w:j + w];

                    // Calculate the mean and standard deviation of the neighbourhood region
                    float wBmn = (float)Median(neighbours);
                    float wBstd = (float)standardDeviation(neighbours);

                    // Calculate the threshold value
                    float wBTH = wBmn * (1 - k * (1-wBstd/R));

                    // Threshold the pixel
                    float aaa = bitmapDataIn[i + j * data.Stride] +
                        bitmapDataIn[i + j * data.Stride + 1] +
                        bitmapDataIn[i + j * data.Stride + 2];
                    aaa /= 3;
                    bitmapDataout[i + j * data.Stride] =
                        bitmapDataout[i + j * data.Stride + 1] =
                        bitmapDataout[i + j * data.Stride + 2] =
                        aaa < wBTH ? byte.MinValue : byte.MaxValue;
                }
            }

            Marshal.Copy(bitmapDataout, 0, data.Scan0, bitmapDataout.Length);
            bitmap.UnlockBits(data);

            return bitmap;
        }
        public Bitmap Phansalkar(Bitmap bitmap, int w = 2, float k = 0.2f, int R = 125, float p = 5, int q =10, float e = 2.718f)
        {
            System.Drawing.Imaging.BitmapData data = null;
            byte[] bitmapDataIn = LockBitmap(bitmap, ref data);
            byte[] bitmapDataout = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapDataout, 0, bitmapDataout.Length);

            int dy = data.Height, dx = data.Stride;
            //int imgSize = data.Height * data.Width;
            //imgN = copy(img);

            // Calculate the radius of the neighbourhood
            //int w = (n - 1) / 2;

            // Process the image
            for (int i = w + 1; i < dx - w; i++)
            {
                for (int j = w + 1; j < dy - w; j++)
                {
                    List<double> neighbours = new List<double>();
                    // Extract the neighbourhood area
                    for (int x = i - w; x < i + w; x++)
                    {
                        for (int y = j - w; y < j + w; y++)
                        {
                            float bbb = bitmapDataIn[x + y * data.Stride] + bitmapDataIn[x + y * data.Stride + 1] + bitmapDataIn[x + y * data.Stride + 2];
                            bbb /= 3;
                            neighbours.Add(bbb);
                        }
                    }
                    //block = bitmapData[i - w:i + w, j - w:j + w];

                    // Calculate the mean and standard deviation of the neighbourhood region
                    float wBmn = (float)Median(neighbours);
                    float wBstd = (float)standardDeviation(neighbours);

                    // Calculate the threshold value
                    float wBTH = (float)(wBmn * ( 1 + p * Math.Pow(e,-q*wBmn) + k * (wBstd/R - 1)));

                    // Threshold the pixel
                    float aaa = bitmapDataIn[i + j * data.Stride] +
                        bitmapDataIn[i + j * data.Stride + 1] +
                        bitmapDataIn[i + j * data.Stride + 2];
                    aaa /= 3;
                    bitmapDataout[i + j * data.Stride] =
                        bitmapDataout[i + j * data.Stride + 1] =
                        bitmapDataout[i + j * data.Stride + 2] =
                        aaa < wBTH ? byte.MinValue : byte.MaxValue;
                }
            }

            Marshal.Copy(bitmapDataout, 0, data.Scan0, bitmapDataout.Length);
            bitmap.UnlockBits(data);

            return bitmap;
        }

        public Bitmap Bernsen(Bitmap bitmap, int w = 2, int l = 15)
        {
            System.Drawing.Imaging.BitmapData data = null;
            byte[] bitmapDataIn = LockBitmap(bitmap, ref data);
            byte[] bitmapDataout = new byte[data.Stride * data.Height];

            Marshal.Copy(data.Scan0, bitmapDataout, 0, bitmapDataout.Length);
            for (int i = 0; i < bitmapDataout.Length; i++)
            {
                bitmapDataout[i] = byte.MaxValue;
            }
            int dy = data.Height, dx = data.Stride;
            //int imgSize = data.Height * data.Width;
            //imgN = copy(img);

            // Calculate the radius of the neighbourhood
            //int w = (n - 1) / 2;

            // Process the image
            for (int i = w + 1; i < dx - w; i++)
            {
                for (int j = w + 1; j < dy - w; j++)
                {
                    List<double> neighbours = new List<double>();
                    // Extract the neighbourhood area
                    for (int x = i - w; x < i + w; x++)
                    {
                        for (int y = j - w; y < j + w; y++)
                        {
                            float bbb = bitmapDataIn[x + y * data.Stride] + bitmapDataIn[x + y * data.Stride + 1] + bitmapDataIn[x + y * data.Stride + 2];
                            bbb /= 3;
                            neighbours.Add(bbb);
                        }
                    }
                    //block = bitmapData[i - w:i + w, j - w:j + w];

                    // Calculate the mean and standard deviation of the neighbourhood region
                    //float wBmn = (float)Median(neighbours);
                    //float wBstd = (float)standardDeviation(neighbours);
                    float wMin = (float)neighbours.Min();
                    float wMax = (float)neighbours.Max();

                    // Calculate the threshold value
                    //float wBTH = (wBmn + k * wBstd);
                    float wBTH = (wMin + wMax) / 2;
                    float localContrast = wMax - wMin;

                    int index = i + j * data.Stride;

                    // Threshold the pixel
                    float aaa = bitmapDataIn[index] +
                        bitmapDataIn[index + 1] +
                        bitmapDataIn[index + 2];
                    aaa /= 3;

                   // // A
                   // bitmapDataout[index] =
                   //     bitmapDataout[index + 1] =
                   //     bitmapDataout[index + 2] =
                   //     aaa < wBTH ? byte.MinValue : byte.MaxValue;
                   //
                   // B
                   if (localContrast < l)
                   {
                       bitmapDataout[index] =
                       bitmapDataout[index + 1] =
                       bitmapDataout[index + 2] =
                       wBTH >= 128 ? byte.MaxValue : byte.MinValue;
                   }
                   else
                   {
                       bitmapDataout[index] =
                       bitmapDataout[index + 1] =
                       bitmapDataout[index + 2] =
                       aaa > wBTH ? byte.MaxValue : byte.MinValue;
                   }

                    // Threshold the pixel

                }
            }

            Marshal.Copy(bitmapDataout, 0, data.Scan0, bitmapDataout.Length);

            bitmap.UnlockBits(data);

            return bitmap;
        }
    }
}
