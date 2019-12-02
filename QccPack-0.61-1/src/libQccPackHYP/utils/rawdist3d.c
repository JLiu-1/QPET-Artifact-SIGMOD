/*
 *
 * QccPack: Quantization, compression, and coding utilities
 * Copyright (C) 1997-2016  James E. Fowler
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */


#include "rawdist3d.h"


#define USG_STRING "[-vo %:] [-mse %:] [-snr %:] [-mae %:] [-s %:] [-bpv %d:bits_per_voxel] [-bsq %:] [-bil %:] [-bip %:] [-bend %:] [-lend %:] %d:num_rows %d:num_cols %d:num_frames %s:rawfile1 %s:rawfile2"

#define ALL 0
#define MSE 1
#define SNR 2
#define MAE 3

int value_only_flag = 0;
int mse_specified = 0;
int snr_specified = 0;
int mae_specified = 0;

int DistortionMethod = ALL;

QccString RawFile1;
QccString RawFile2;

int NumRows;
int NumCols;
int NumFrames;

int Signed = 0;

int BitsPerVoxel = 16;

int Format;
int FormatBSQSpecified = 0;
int FormatBILSpecified = 0;
int FormatBIPSpecified = 0;

int Endian;
int BigEndianSpecified = 0;
int LittleEndianSpecified = 0;


int main(int argc, char *argv[])
{
  double value = 0;
  double mse;
  double mae;
  double snr;

  QccInit(argc, argv);
  
  if (QccParseParameters(argc, argv,
			 USG_STRING,
                         &value_only_flag,
                         &mse_specified,
                         &snr_specified,
                         &mae_specified,
                         &Signed,
                         &BitsPerVoxel,
                         &FormatBSQSpecified,
                         &FormatBILSpecified,
                         &FormatBIPSpecified,
                         &BigEndianSpecified,
                         &LittleEndianSpecified,
                         &NumRows,
                         &NumCols,
                         &NumFrames,
                         RawFile1,
                         RawFile2))
    QccErrorExit();
  
  if (mse_specified + snr_specified + mae_specified > 1)
    {
      QccErrorAddMessage("%s: Only one of -mse, -snr, or -mae can be specified",
                         argv[0]);
      QccErrorExit();
    }

  if (mse_specified)
    DistortionMethod = MSE;
  else if (snr_specified)
    DistortionMethod = SNR;
  else if (mae_specified)
    DistortionMethod = MAE;
  else
    DistortionMethod = ALL;

  if (value_only_flag &&
      (!mse_specified && !snr_specified && !mae_specified))
    {
      QccErrorAddMessage("%s: One of -mse, -snr, and -mae must be used with -vo",
                         argv[0]);
      QccErrorExit();
    }

  if (FormatBSQSpecified + FormatBILSpecified + FormatBIPSpecified > 1)
    {
      QccErrorAddMessage("%s: Only one of -bsq, -bil, or -bip can be specified",
                         argv[0]);
      QccErrorExit();
    }

  if (FormatBSQSpecified)
    Format = QCCHYP_RAWFORMAT_BSQ;
  else
    if (FormatBILSpecified)
      Format = QCCHYP_RAWFORMAT_BIL;
    else
      if (FormatBIPSpecified)
        Format = QCCHYP_RAWFORMAT_BIP;
      else
        Format = QCCHYP_RAWFORMAT_BSQ;

  if ((BitsPerVoxel < 1) ||
      (BitsPerVoxel > 32))
    {
      QccErrorAddMessage("%s: bpv must be between 1 and 32",
                         argv[0]);
      QccErrorExit();
    }
      
  if (BigEndianSpecified + LittleEndianSpecified > 1)
    {
      QccErrorAddMessage("%s: Only one of -be or -le can be specified",
                         argv[0]);
      QccErrorExit();
    }

  if (BigEndianSpecified)
    Endian = QCCHYP_RAWENDIAN_BIG;
  else
    if (LittleEndianSpecified)
      Endian = QCCHYP_RAWENDIAN_LITTLE;
    else
      Endian = QCCHYP_RAWENDIAN_BIG;

  if (QccHYPRawDist3D(RawFile1, RawFile2,
                      &mse, &mae, &snr,
                      NumFrames, NumRows, NumCols,
                      BitsPerVoxel, Signed, Format, Endian))
    {
      QccErrorAddMessage("%s: Error calling QccHYPRawDist3D()",
                         argv[0]);
      QccErrorExit();
    }

  if (!value_only_flag)
    printf("The distortion between files\n   %s\n         and\n   %s\nis:\n",
           RawFile1, RawFile2);
  
  if (!value_only_flag)
    {
      switch (DistortionMethod)
        {
        case MSE:
          printf("%16.6f (MSE)\n", mse);
          break;
        case SNR:
          printf("%16.6f dB (SNR)\n", snr);
          break;
        case MAE:
          printf("%16.6f (MAE)\n", mae);
          break;
        case ALL:
        default:
          printf("%16.6f\t(MSE)\n%16.6f dB\t(SNR)\n%16.6f\t(MAE)\n",
                 mse, snr, mae);
          break;
        }
    }
  else
    {
      switch (DistortionMethod)
        {
        case MSE:
          value = mse;
          break;
        case SNR:
          value = snr;
          break;
        case MAE:
          value = mae;
          break;
        }
      printf("%f\n", value);
    }

  QccExit;
}
