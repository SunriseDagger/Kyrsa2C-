using System;
using System.IO;
using System.Text;
using System.Linq;
using System.Numerics;

namespace Kyrsa2
{
	class Program
	{
		static double EPS = 1.0e-14; //погрешность
		static double OBS = 1;//1.0e+16; //1; //для увеличения диагонали и краевых условий

		// количество x и y
		static string[] ss  = File.ReadAllLines("C:\\Users\\МSI\\source\\repos\\Kyrsa2\\Kyrsa2\\info.txt");
		static int Nx = int.Parse(ss[0]);
		static int Ny = int.Parse(ss[1]);
		static int Nt = int.Parse(ss[2]); //количество временных слоёв (без 0)
		static int nS1 = int.Parse(ss[3]); //список узлов, лежащих на первых краевых
		static double[] AllX = { 0, 1, 2}; //сетка по х
		static double[] AllY = { 0, 1, 2}; //сетка по y

		//static double[] AllT = { 0, 0.1, 0.3, 0.7, 1 }; //сетка по времени
		//static double[] AllT = { 0, 1, 2, 3 }; //сетка по времени
		static double[] AllT = { 0, 0.5, 1, 1.5, 2 }; //сетка по времени
		//static double[] AllT = { 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2};
		//static double[] AllT = { 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3 };
		static int N = Nx * Ny; //размерность матрицы, количество узлов
		static int feN = (Nx - 1) * (Ny - 1); //количество конечных элементов

		static string[] sS1 = File.ReadAllLines("C:\\Users\\МSI\\source\\repos\\Kyrsa2\\Kyrsa2\\S1.txt");
		static int[] S1 = new int[sS1.Length];

		static double[] di = new double[N]; //диагональ глобальной матрицы 
		static int[] ig; //количество элементов гм
		static int[] jg; //позиции элементов гм

		static double[] ggl; //нижний треугольник гм
		static double[] ggu; //верхний треугольник гм

		static double[] Gdi = new double[N]; //диагональ ж
		static double[] Gggl; //нижний треугольник ж
		static double[] Gggu; //верхний треугольник ж
		static double[] Mdi = new double[N]; //диагональ м
		static double[] Mggl; //нижний треугольник м
		static double[] Mggu; //верхний треугольник м

		static double[] fv = new double[4];     // значения правой части
		static double[,] G = new double[4, 4];   // матрица жёсткости
		static double[,] M = new double[4, 4];   // матрица масс
		static double[] F = new double[4];      // вектор правой части
		static int[] Li = new int[4];      // для занесения в глобальную

		static double[] L;  // нижний треугольник обусловленной матрицы
		static double[] U; // верхний треугольник обусловленной матрицы
		static double[] diag = new double[N]; // диагональ обусловленной матрицы
		static double[] x = new double[N]; // вектор решения
		static double[] b1 = new double[N]; //вектор правой части
		static double[] b0 = new double[N]; //вектор предыдущей правой части
		static double[] q0 = new double[N]; //вектор предыдущего решения
		static double[] d = new double[N]; // общий вектор правой части		

		static int[,] gg = new int[2, N];

		static void Main(string[] args)
		{
			for (int i = 0; i < sS1.Length; i++)
            {
				S1[i] = int.Parse(sS1[i]);
            }

			int yy = 0;
			int z = 0;

			for (int xx = 0; xx < Nx; xx++)
            {
				gg[0, z] = xx; //x
				gg[1, z] = yy; //y
				z++;
				if (xx + 1 == Nx)
                {
					xx = -1;
					yy++;
                }

				if (z == N) break;
			}

			int n = 0;

			for (int k = 0; k < Ny; k++)
            {
				for (int i = 0; i < Nx; i++)
				{
					q0[n] = Uf(AllX[i], AllY[k], 0);
					b0[n] = Ff(AllX[i], AllY[k], 0);
					n++;
				}
			}

			GeneratePortret();
			main1();
			Console.ReadLine();
		}

		static void main1()
		{
			double dt;
			for (int i = 1; i < Nt; i++)
            {
				dt = AllT[i] - AllT[i - 1];
				localMatrix(AllT[i]);
				assemblyA(dt);

				double[] y = new double[N];
				MultMonV(q0, y, dt);
				
				for (int j = 0; j < N; j++) d[j] = ((b1[j] + b0[j]) / 2.0) + y[j];
				b1 = d;

				//делаю первые краевые
				for (int k = 0; k < N; k++)
                {
					if (S1.Contains(k)) Bound1(k, AllT[i]);
				}
				
				Console.WriteLine("t = " + AllT[i]);
				run();
				q0 = result(AllT[i]);
				b0 = b1;
			}
			
		}

		static double Ff(double x, double y, double t) //f
		{
			return 2*t;
			//return x + y;
			//return 75 * Math.Pow(t, 14) + 2;
			//return 15 * Math.Pow(t, 2) * x + 2;
			//return 10 * t * x + 2;
			//return 5 * x + y;
			//return Math.Pow(x, 2) + Math.Pow(y, 2);
			//return 4 * x + 4 * y;
			//return 3 * Math.Pow(t, 2);
		}

		static double Uf(double x, double y, double t) //u* - точное
		{
			//return 1;
			//return x + y;
			//return Math.Pow(x, 4);
			//return 5 * Math.Pow(t, 2) * x + 2 * t;
			//return 5 * x * t + 2 * t;
			//return 5 * Math.Pow(x, 3) * Math.Pow(t, 2) + 2 * t;
			//return t;
			return Math.Pow(t, 2);
		}

		static void GeneratePortret()
        {
			int[] listbeg = new int[N];
			for (int i = 0; i < N; i++) listbeg[i] = -1;

			int[] list1 = new int[N * N];
			int[] list2 = new int[N * N];
			int listsize = -1, iaddr, ind1, ind2, k;

			for (int iel = 0; iel < feN; iel++) //проход по всем кэ
            {
				for (int i = 0; i < 4; i++) //по 4м узлам кэ
                {
					k = knots_num(iel, i); //находим глобальный номер узла
					for (int j = i + 1; j < 4; j++)
                    {
						ind1 = k;
						ind2 = knots_num(iel, j);

						if (ind2 < ind1) //занесение связи большого номера с меньшим, те ind2 c ind1
                        {
							ind1 = ind2;
							ind2 = k;
                        }
						iaddr = listbeg[ind2];

						if (iaddr == -1) //если список пуст
                        {
							listsize++;
							listbeg[ind2] = listsize;
							list1[listsize] = ind1;
							list2[listsize] = -1;
                        }
						else //ищем в списке ind1
                        {
							while (list1[iaddr] < ind1 && list2[iaddr] > -1) iaddr = list2[iaddr];
							if (list1[iaddr] > ind1) //если не нашли и встретили эл с большим номером
                            {						 //добавляем перед ним, чтобы список был упорядоченным
								listsize++;
								list1[listsize] = list1[iaddr];
								list2[listsize] = list2[iaddr];
								list1[iaddr] = ind1;
								list2[iaddr] = listsize;
							}
							else if (list1[iaddr] < ind1) //если не нашли, то вставляем в конец
							{                       
								listsize++;
								list2[iaddr] = listsize;
								list1[listsize] = ind1;
								list2[listsize] = -1;
							}
						}
					}
                }
            }

			ig = new int[N + 1];
			jg = new int[listsize + 1];
			int nt = jg.Length;

			ggl = new double[nt]; //нижний треугольник гм
			ggu = new double[nt]; //верхний треугольник гм

			Gggl = new double[nt]; //нижний треугольник ж
			Gggu = new double[nt]; //верхний треугольник ж

			Mggl = new double[nt]; //нижний треугольник м
			Mggu = new double[nt]; //верхний треугольник м

			L = new double[nt];  // нижний треугольник обусловленной матрицы
			U = new double[nt]; // верхний треугольник обусловленной матрицы

			for (int i = 0; i < N; i++)
            {
				ig[i + 1] = ig[i];

				for (iaddr = listbeg[i]; iaddr != -1; )
                {
					jg[ig[i + 1]] = list1[iaddr];
					ig[i + 1]++;
					iaddr = list2[iaddr];
                }
            }

		}

		static int knots_num(int el, int localKnots)
        {
			int globalKnonts = 0, k = 0, i = 0;

			switch (el)
            {
				case 0:
					k = 0;
					i = 0;
					break;
				case 1:
					k = 0;
					i = 1;
					break;
				case 2:
					k = 1;
					i = 0;
					break;
				case 3:
					k = 1;
					i = 1;
					break;
				/*case 4:
					k = 1;
					i = 1;
					break;
				case 5:
					k = 1;
					i = 2;
					break;
				case 6:
					k = 2;
					i = 0;
					break;
				case 7:
					k = 2;
					i = 1;
					break;
				case 8:
					k = 2;
					i = 2;
					break;*/
				/*case 9:
					k = 2;
					i = 1;
					break;
				case 10:
					k = 2;
					i = 2;
					break;
				case 11:
					k = 2;
					i = 3;
					break;
				case 12:
					k = 3;
					i = 0;
					break;
				case 13:
					k = 3;
					i = 1;
					break;
				case 14:
					k = 3;
					i = 2;
					break;
				case 15:
					k = 3;
					i = 3;
					break;*/
			}

			switch (localKnots)
            {
				case 0:
					globalKnonts = Nx * k + i;
					break;

				case 1:
					globalKnonts = Nx * k + i + 1;
					break;

				case 2:
					globalKnonts = Nx * (k + 1) + i;
					break;

				case 3:
					globalKnonts = Nx * (k + 1) + i + 1;
					break;
			}

			return globalKnonts;
        }

		static void Bound1(int i, double t)
		{
			int xx, yy;
			xx = gg[0, i];
			yy = gg[1, i];

			int k, j;
			di[i] = 1;
			b1[i] = Uf(AllX[xx], AllY[yy], t);

			for (k = ig[i]; k < ig[i + 1]; k++) ggl[k] = 0;

			for (k = i + 1; k < N; k++)
            {
				for (j = ig[k]; j < ig[k + 1]; j++)
                {
					if (jg[j] == i) ggu[j] = 0;
                }
			}
		}

		static double[] MultMonV (double[] a, double[] y, double dt) //для правой части, где матрицу умножаю на q^(i-1)
        {
			double dt1 = 1 / dt;
			int i, j;
			double[] b = new double[N];

			for (i = 0; i < N; i++)
			{
				b[i] = (dt1 * Mdi[i] - 0.5 * Gdi[i]) * a[i];
				for (j = ig[i]; j < ig[i + 1]; j++)
				{
					int z = jg[j];
					b[i] += (dt1 * Mggl[j] - Gggl[j] / 2.0) * a[z];
					b[z] += (dt1 * Mggu[j] - Gggl[j] / 2.0) * a[i];
				}
			}

			for (i = 0; i < N; i++) y[i] = b[i];

			return y;
		}

		static void localMatrix(double dt) //создание локальной матрицы, для последующего присоединения к глобальной
		{
			double lambda, sigma, px, y, xp, yp, hx, hy;  // параметры КЭ
			double ud, u1, u2, u3;

			// коэффициенты
			string[] s = File.ReadAllLines("C:\\Users\\МSI\\source\\repos\\Kyrsa2\\Kyrsa2\\math.txt");
			lambda = double.Parse(s[0]);
			sigma = double.Parse(s[1]);

			//Создание элементов локальной матрицы
			for (int k = 0; k < Ny - 1; k++)
				for (int i = 0; i < Nx - 1; i++)
				{

					px = AllX[i];  // опорная точка
					y = AllY[k];

					// верхние точки
					xp = AllX[i + 1];
					yp = AllY[k + 1];

					// шаги
					hx = xp - px;
					hy = yp - y;


					// значения правой части
					fv[0] = Ff(px, y, dt);
					fv[1] = Ff(xp, y, dt);
					fv[2] = Ff(px, yp, dt);
					fv[3] = Ff(xp, yp, dt); //dt=t1-t0

					// Задаём значения матрицы жёсткости G
					ud = (hy / hx + hx / hy) / 3.0;
					u1 = (-hy / hx + hx / (2.0 * hy)) / 3.0;
					u2 = (hy / (2.0 * hx) - hx / hy) / 3.0;
					u3 = (hy / hx + hx / hy) / -6f;


					G[0, 0] = G[1, 1] = G[2, 2] = G[3, 3] = ud;
					G[0, 1] = G[1, 0] = G[2, 3] = G[3, 2] = u1;
					G[0, 2] = G[1, 3] = G[2, 0] = G[3, 1] = u2;
					G[0, 3] = G[1, 2] = G[2, 1] = G[3, 0] = u3;

					// Задаём значения матрицы масс M
					ud = hx * hy / 9.0;
					u1 = hx * hy / 18.0;
					u2 = hx * hy / 36.0;

					M[0, 0] = M[1, 1] = M[2, 2] = M[3, 3] = ud;
					M[0, 1] = M[1, 0] = M[0, 2] = M[2, 0] = M[1, 3] = M[3, 1] = M[2, 3] = M[3, 2] = u1;
					M[0, 3] = M[1, 2] = M[2, 1] = M[3, 0] = u2;

					for (int q = 0; q < 4; q++)
                    {
						F[q] = M[q, 0] * fv[0] + M[q, 1] * fv[1] + M[q, 2] * fv[2] + M[q, 3] * fv[3];
					}

					Li[0] = Nx * k + i;
					Li[1] = Nx * k + i + 1;
					Li[2] = Nx * (k + 1) + i;
					Li[3] = Nx * (k + 1) + i + 1;

					//Console.WriteLine("i = " + i + ", k = " + k + ", F[0] = " + F[0] + ", F[1] = " + F[1] + ", F[2] = " + F[2] + ", F[3] = " + F[3]);
					for (int z = 0; z < 4; z++)
					{
						b1[Li[z]] += F[z];
					}
					/*Console.WriteLine("i = " + i + ", k = " + k);
					for (int m = 0; m < 4; m++)
                    {
						Console.WriteLine(m);
						Console.WriteLine("G[0,0] = " + G[m, 0] + "G[0,1] = " + G[m, 1] + "G[0,2] = " + G[m, 2] + "G[0,3] = " + G[m, 3]);
						Console.WriteLine("M[0,0] = " + M[m, 0] + "M[0,1] = " + M[m, 1] + "M[0,2] = " + M[m, 2] + "M[0,3] = " + M[m, 3]);
					}*/
					

					PortrG(G);
					PortrM(M);

				}

			//умножаем жесткость на лямбду и массу на сигму
			for (int i = 0; i < di.Length; i++)
            {
				Gdi[i] *= lambda;
				Mdi[i] *= sigma;
			}

			for (int i = 0; i < ggu.Length; i++)
			{
				Gggu[i] *= lambda;
				Mggu[i] *= sigma;
				Gggl[i] *= lambda;
				Mggl[i] *= sigma;
			}	

		}

		//	Добавление элементов в портрет матрицы
		static void PortrG (double[,] A)
		{
			int ibeg, iend, ind;
			//double[,] A = new double[4, 4];
			//сюда отправляется L(k), A(k,k)
			//заносим диагональные элементы
			for (int i = 0; i < 4; i++)
            {
				Gdi[Li[i]] += A[i, i];
            }

			for (int i = 0; i < 4; i++)
			{
				ibeg = ig[Li[i]];

				for (int j = 0; j < i; j++)
				{
					iend = ig[Li[i] + 1] - 1;

					while (jg[ibeg] != Li[j])
					{
						ind = (ibeg + iend) / 2 + 1;

						if (jg[ind] <= Li[j])
						{
							ibeg = ind;
						}
						else
						{
							//iend = ind;
							iend = (iend == ind) ? ind - 1 : ind;
						}
					}

					Gggl[ibeg] += A[i, j];
					Gggu[ibeg] += A[j, i];
					ibeg++;
				}
			}

		}
		static void PortrM (double[,] A)
		{
			int ibeg, iend, ind;
			//double[,] A = new double[4, 4];
			//сюда отправляется L(k), A(k,k)
			//заносим диагональные элементы
			for (int i = 0; i < 4; i++)
			{
				Mdi[Li[i]] += A[i, i];
			}

			for (int i = 0; i < 4; i++)
			{
				ibeg = ig[Li[i]];

				for (int j = 0; j < i; j++)
				{
					iend = ig[Li[i] + 1] - 1; 

					while (jg[ibeg] != Li[j])
					{
						/*Console.WriteLine("iend =" + iend);
						Console.WriteLine("ibeg =" + ibeg);*/
						ind = (ibeg + iend) / 2 + 1;
						//Console.WriteLine("ind =" + ind);
						if (jg[ind] <= Li[j])
						{
							ibeg = ind;
						}
						else
						{
							//iend = ind;
							iend = (iend == ind) ? ind - 1 : ind;
						}
					}

					Mggl[ibeg] += A[i, j];
					Mggu[ibeg] += A[j, i];
					ibeg++;
				}
			}
		}

		static void assemblyA(double dt)
        {
			double dt1 = 1 / dt;

			for (int i = 0; i < N; i++)
            {
				di[i] = dt1 * Mdi[i] + Gdi[i] / 2;
			
			}

			for (int i = 0; i < ggu.Length; i++)
            {
				ggu[i] = dt1 * Mggu[i] + Gggu[i] / 2.0;
				ggl[i] = dt1 * Mggl[i] + Gggl[i] / 2.0;
			}
        }

		//	Разложение LU
		static int LU()
		{

			int i, j;
			for (i = 0; i < N; i++)
			{
				for (j = ig[i]; (j < ig[i + 1]); j++)
				{
					L[j] = (ggl[j] - MultLU(i, jg[j]));
					U[j] = (ggu[j] - MultLU(jg[j], i)) / diag[jg[j]];
				}
				diag[i] = di[i] - MultLU(i, i);
			}
			return 0;
		}

		//	Решение нижнего треугольника, приближение
		static void solL(double[] a, double[] b)
		{
			int i, j;
			double result;

			for (i = 0; i < N; i++)
			{
				result = 0;
				for (j = ig[i]; (j < ig[i + 1]); j++)
				{
					int z = jg[j];
					result += L[j] * b[z];
				}
				b[i] = (a[i] - result) / diag[i];
			}
		}

		//	Решение верхнего треугольника
		static void solU(double[] a, double[] b)
		{
			int i, j;
			for (i = 0; i < N; i++) b[i] = a[i];
			for (i = N - 1; i >= 0; i--)
			{
				for (j = ig[i]; (j < ig[i + 1]); j++)
				{
					int z = jg[j];
					b[z] -= U[j] * b[i];
				}
			}
		}

		//	Умножение матрицы на вектор
		static void MultMatrixOnVector(double[] a, double[] b)
		{
			int i, j;
			double[] bl = new double[N];
			for (i = 0; i < N; i++)
			{
				bl[i] = di[i] * a[i];
				for (j = ig[i]; (j < ig[i + 1]); j++)
				{
					int z = jg[j];
					bl[i] += ggl[j] * a[z];
					bl[z] += ggu[j] * a[i];
				}
			}

			for (i = 0; i < N; i++)
				b[i] = bl[i];
		}


		//	Скалярное произведение двух векторов
		static double ScalarMult(double[] a, double[] b)
		{
			int i;
			double result = 0;
			for (i = 0; i < N; i++)
			{
				result += a[i] * b[i];
			}
			return result;
		}

		//	Скалярное произведение для факторизации
		static double MultLU(int i, int j)
		{
			int k, l, find;
			double result = 0.0;
			if (i == j)
			{
				for (k = ig[i]; (k < ig[i + 1]); k++)
					result += U[k] * L[k];
			}
			else
			{
				// верхний треугольник
				if (i > j)
				{
					for (k = ig[j]; k < ig[j + 1]; k++)
					{
						find = 0;
						for (l = ig[i]; l < (ig[i + 1]) && find == 0; l++)
						{
							if (jg[l] == jg[k])
							{
								result += U[k] * L[l];
								find = 1;
							}
						}
					}
				}
				// нижний треугольник
				else
				{
					for (l = ig[i]; (l < ig[i + 1]); l++)
					{
						find = 0;
						for (k = ig[j]; k < ig[j + 1] && find == 0; k++)
						{
							if (jg[l] == jg[k])
							{
								result += U[k] * L[l];
								find = 1;
							}
						}
					}
				}

			}
			return result;
		}

		static double NormV (double[] a) //норма вектора
        {
			double n = 0;
			for (int i = 0; i < a.Length; i++)
            {
				n += Math.Pow(a[i], 2);
            }

			Math.Pow(n, 0.5);

			return n;
        }

		//	Основная функция метода
		static void run()
		{
			int iter;
			int check, stop;
			double alpha, alphazn, alphach, beta, betach, betazn, CheckExit;

			double[] r0 = new double[N]; // вспомогательный вектор
			double[] r = new double[N]; // вспомогательный вектор
			double[] p = new double[N]; // вспомогательный вектор
			double[] z = new double[N]; // вспомогательный вектор
			double[] q = new double[N]; // вспомогательный вектор	
			double[] s = new double[N]; // вспомогательный вектор				
			double[] sout = new double[N]; // вспомогательный вектор

			// Факторизация
			check = LU();
			//если на диагонали нули
			if (check != 0) Console.WriteLine("Нельзя выполнить факторизацию", check + 1); 

			// Инициализация
			stop = 0;
			for (int i = 0; i < N; i++) x[i] = 0; //начальное приближение 
			solL(b1, r);
			solU(r, z);
			MultMatrixOnVector(z, q);
			solL(q, p);

			// Процесс решения СЛАУ
			for (iter = 0; iter < 10000 && stop == 0; iter++)
			{
				alphach = ScalarMult(p, r);
				alphazn = ScalarMult(p, p);
				alpha = alphach / alphazn;
				for (int i = 0; i < N; i++) x[i] += alpha * z[i];
				for (int i = 0; i < N; i++) r[i] -= alpha * p[i];
				solU(r, s);
				MultMatrixOnVector(s, sout);
				solL(sout, q);
				betazn = ScalarMult(p, p);
				betach = ScalarMult(p, q);
				beta = -betach / betazn;
				//solU(r, s);
				for (int i  = 0; i < N; i++) z[i] = beta * z[i] + s[i];
				for (int i  = 0; i < N; i++) p[i] = beta * p[i] + q[i];
				CheckExit = ScalarMult(r, r);
				//CheckExit = NormV(r)/NormV(b1);
				if (CheckExit < EPS) stop = 1;
			}

			Console.WriteLine("iter = " + iter);
		}

		static double[] result(double t)
		{

			int i, k, num = 0;
			double px, y, func, res = 0.0, norm = 0.0, tmp;

			Console.WriteLine("X" + "            " + "Y" + "            " + "U" + "            " + "U*" + "            " + "U*-U");
		

			for (k = 0; k < Ny; k++)
            {
				for (i = 0; i < Nx; i++)
				{
					px = AllX[i];
					y = AllY[k];
					func = Uf(px, y, t);
					tmp = Math.Abs(x[num] - func);
					res += tmp * tmp;
					norm += func * func;
					Console.WriteLine(px + "            " + y + "            " + x[num] + "            " + func + "            " + tmp);
					num++;
				}
			}
				

			Console.WriteLine("\n ||U-U*|| / ||U*|| = " + Math.Sqrt(res / norm));
			return x;
		}
	}
}
