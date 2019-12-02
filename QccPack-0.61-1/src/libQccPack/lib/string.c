/*
 * 
 * QccPack: Quantization, compression, and coding libraries
 * Copyright (C) 1997-2016  James E. Fowler
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 * 
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
 * MA 02139, USA.
 * 
 */


#include "libQccPack.h"


void QccStringMakeNull(QccString qccstring)
{
  if (qccstring == NULL)
    return;

  qccstring[0] = '\0';
}


int QccStringNull(const QccString qccstring)
{
  if (qccstring == NULL)
    return(1);
  return(strlen(qccstring) <= 0);
}


void QccConvertToQccString(QccString qccstring, const char *str)
{
  if ((qccstring == NULL) || (str == NULL))
    return;

  strncpy((char *)qccstring, str, QCCSTRINGLEN);
  qccstring[QCCSTRINGLEN] = '\0';
}


void QccStringCopy(QccString qccstring1, const QccString qccstring2)
{
  if (qccstring1 == NULL)
    return;
  if (qccstring2 == NULL)
    return;

  strncpy((char *)qccstring1, (char *)qccstring2, QCCSTRINGLEN);
  qccstring1[QCCSTRINGLEN] = '\0';
  
}


void QccStringSprintf(QccString qccstring, const char *format, ...)
{
  va_list ap;

  va_start(ap, format);

#ifdef QCC_NO_SNPRINTF
  vsprintf(qccstring, format, ap);
#else
  vsnprintf(qccstring, QCCSTRINGLEN, format, ap);
#endif

  qccstring[QCCSTRINGLEN] = '\0';

  va_end(ap);
}
