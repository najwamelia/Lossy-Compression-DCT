# -*- coding:utf-8 -*-
from PIL import Image
import numpy as np
from ctypes import c_ubyte
import math
import sys
import imageio

img = Image.open('original.gif')

dct_blockList = []

# Buat blok 8x8 -> 1024
for i in range(32):
    for j in range(32):
        box = (j * 8, i * 8, (j + 1) * 8, (i + 1) * 8)
        block = img.crop(box)
        blocks = np.array(block)
        dct_blockList.append(blocks)

dct_blockLists = np.zeros((1024, 8, 8))

# Hitung dct untuk setiap blok
for l in range(1024):
    for u in range(8):
        for v in range(8):
            sum = 0.0
            for i in range(8):
                for j in range(8):
                    sum += dct_blockList[l][i][j]*np.cos(
                        ((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
            dct_blockLists[l][u][v] = (1/4)*(Cu*Cv*sum)
            dct_blockLists = dct_blockLists.astype(int)

            f1 = open("dct_result.txt", 'a')
            sys.stdout = f1
            print(dct_blockLists[l])
sys.stdout = sys.__stdout__
f1.close()

# Jalankan kuantisasi untuk setiap blok
quantization_table = [[16, 11, 10, 16, 24, 40, 51, 61],
                      [12, 12, 14, 19, 26, 58, 60, 55],
                      [14, 13, 16, 24, 40, 57, 69, 56],
                      [14, 17, 22, 29, 51, 87, 80, 62],
                      [18, 22, 37, 56, 68, 109, 103, 77],
                      [24, 35, 55, 64, 81, 104, 113, 92],
                      [49, 64, 78, 87, 103, 121, 120, 101],
                      [72, 92, 95, 98, 112, 100, 103, 99]]

quantization_blockList = np.zeros((1024, 8, 8))

for i in range(1024):
    for j in range(8):
        for k in range(8):
            quantization_blockList[i][j][k] = np.around(
                dct_blockLists[i][j][k] / quantization_table[j][k])

    f2 = open("quantized_block_result.txt", 'a')
    sys.stdout = f2
    print(quantization_blockList[i])
sys.stdout = sys.__stdout__
f2.close()

# Jalankan koefisien DCT terkuantisasi terbalik untuk setiap blok
inverse_qunatized_dct_blockList = np.zeros((1024, 8, 8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            inverse_qunatized_dct_blockList[i][j][k] = quantization_blockList[i][j][k] * \
                quantization_table[j][k]

    f3 = open("inverse_quantized_dct_block_result.txt", 'a')
    sys.stdout = f3
    print(inverse_qunatized_dct_blockList[i])
sys.stdout = sys.__stdout__
f3.close()


# idct untuk 3-(a)
a_idct = np.zeros((1024, 8, 8))
for l in range(1024):
    for i in range(8):
        for j in range(8):
            sum = 0
            for u in range(8):
                for v in range(8):
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
                    sum += inverse_qunatized_dct_blockList[l][u][v]*np.cos(
                        ((2*i+1)*u*math.pi)/16)*np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            a_idct[l][i][j] = sum
            a_idct = np.clip(a_idct, 0, 255)
            a_idct = a_idct.astype(c_ubyte)

# merekonstruksi gambar untuk 3-(a)
a_reconstruct_image = Image.new('L', (256, 256))
n = 0
for i in range(0, 256, 8):
    for j in range(0, 256, 8):
        b_image = Image.fromarray(a_idct[n])
        area = (j, i, j+8, i+8)
        a_reconstruct_image.paste(b_image, area)
        n += 1
a_reconstruct_image.save("3(a)_reconstruct_image.jpg")

# 3-(b) pertahankan F(0,0) potong yang lainnya
b = np.zeros((1024, 8, 8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if j is 0 and k is 0:
                b[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                b[i][j][k] = 0

    f4 = open("truncate_a.txt", 'a')
    sys.stdout = f4
    print([i])
sys.stdout = sys.__stdout__
f4.close()

# 3-(c) pertahankan F(0,0), F(0,1), F(1,0) potong semua yang lain
c = np.zeros((1024, 8, 8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if (j is 0 and k is 0) or (j is 0 and k is 1) or (j is 1 and k is 0):
                c[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                c[i][j][k] = 0

    f5 = open("truncate_b.txt", 'a')
    sys.stdout = f5
    print(b[i])
sys.stdout = sys.__stdout__
f5.close()


# 3-(d) pertahankan F(0,0), F(1,0), F(0,1), F(1,1), F(2,0), F(0,2) semuanya memotong
d = np.zeros((1024, 8, 8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if (j is 0 and k is 0) or (j is 1 and k is 0) or \
                    (j is 0 and k is 1) or (j is 1 and k is 1) or (j is 2 and k is 0) or (j is 0 and k is 2):
                d[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                d[i][j][k] = 0

    f6 = open("truncate_c.txt", 'a')
    sys.stdout = f6
    print(c[i])
sys.stdout = sys.__stdout__
f6.close()


# 3-(d) F(0,0), F(1,0), F(0,1), F(1,1), F(2,0), F(0,2), F( Jauhkan 0,3), F(1,2), F(2,1), F(3,0) potong semua yang lain
e = np.zeros((1024, 8, 8))
for i in range(1024):
    for j in range(8):
        for k in range(8):
            if (j is 0 and k is 0) or (j is 1 and k is 0) or\
                    (j is 0 and k is 1) or (j is 1 and k is 1) or (j is 2 and k is 0) or\
                    (j is 0 and k is 0) or (j is 1 and k is 2) or (j is 2 and k is 1) or\
                    (j is 3 and k is 0) or (j is 0 and k is 3):
                e[i][j][k] = inverse_qunatized_dct_blockList[i][j][k]
            else:
                e[i][j][k] = 0

    f7 = open("truncate_d.txt", 'a')
    sys.stdout = f7
    print(d[i])
sys.stdout = sys.__stdout__
f7.close()


# 3 - idct untuk (b)
b_idct = np.zeros((1024, 8, 8))
for l in range(1024):
    for i in range(8):
        for j in range(8):
            sum = 0
            for u in range(8):
                for v in range(8):
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
                    sum += b[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16) * \
                        np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            b_idct[l][i][j] = sum
            b_idct = np.clip(b_idct, 0, 255)
            b_idct = b_idct.astype(c_ubyte)

# merekonstruksi gambar untuk 3-(b)
b_reconstruct_image = Image.new('L', (256, 256))
n = 0
for i in range(0, 256, 8):
    for j in range(0, 256, 8):
        b_image = Image.fromarray(b_idct[n])
        area = (j, i, j+8, i+8)
        b_reconstruct_image.paste(b_image, area)
        n += 1
b_reconstruct_image.save("3(b)_reconstruct_image.jpg")

# idct untuk 3-(c)
c_idct = np.zeros((1024, 8, 8))
for l in range(1024):
    for i in range(8):
        for j in range(8):
            sum = 0
            for u in range(8):
                for v in range(8):
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
                    sum += c[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16) * \
                        np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            c_idct[l][i][j] = sum
            c_idct = np.clip(c_idct, 0, 255)
            c_idct = c_idct.astype(c_ubyte)

# merekonstruksi gambar untuk 3-(c)
c_reconstruct_image = Image.new('L', (256, 256))
n = 0
for i in range(0, 256, 8):
    for j in range(0, 256, 8):
        c_image = Image.fromarray(c_idct[n])
        area = (j, i, j+8, i+8)
        c_reconstruct_image.paste(c_image, area)
        n += 1
c_reconstruct_image.save("3(c)_reconstruct_image.jpg")

# idct untuk 3-(d)
d_idct = np.zeros((1024, 8, 8))
for l in range(1024):
    for i in range(8):
        for j in range(8):
            sum = 0
            for u in range(8):
                for v in range(8):
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
                    sum += d[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16) * \
                        np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            d_idct[l][i][j] = sum
            d_idct = np.clip(d_idct, 0, 255)
            d_idct = d_idct.astype(c_ubyte)

# merekonstruksi gambar untuk 3-(d)
d_reconstruct_image = Image.new('L', (256, 256))
n = 0
for i in range(0, 256, 8):
    for j in range(0, 256, 8):
        d_image = Image.fromarray(d_idct[n])
        area = (j, i, j+8, i+8)
        d_reconstruct_image.paste(d_image, area)
        n += 1
d_reconstruct_image.save("3(d)_reconstruct_image.jpg")

# idct untuk 3-(e)
e_idct = np.zeros((1024, 8, 8))
for l in range(1024):
    for i in range(8):
        for j in range(8):
            sum = 0
            for u in range(8):
                for v in range(8):
                    if u == 0:
                        Cu = 1 / math.sqrt(2)
                    else:
                        Cu = 1
                    if v == 0:
                        Cv = 1 / math.sqrt(2)
                    else:
                        Cv = 1
                    sum += e[l][u][v]*np.cos(((2*i+1)*u*math.pi)/16) * \
                        np.cos(((2*j+1)*v*math.pi)/16)*(1/4)*Cu*Cv
            e_idct[l][i][j] = sum
            e_idct = np.clip(e_idct, 0, 255)
            e_idct = e_idct.astype(c_ubyte)

# merekonstruksi gambar untuk 3-(e)
e_reconstruct_image = Image.new('L', (256, 256))
n = 0
for i in range(0, 256, 8):
    for j in range(0, 256, 8):
        e_image = Image.fromarray(e_idct[n])
        area = (j, i, j+8, i+8)
        e_reconstruct_image.paste(e_image, area)
        n += 1
e_reconstruct_image.save("3(e)_reconstruct_image.jpg")

# 3-(a) Hitung MSE
MSE_A = 0
MSE_A = np.square(np.subtract(dct_blockList, a_idct)).mean(axis=None)
print("3-(a) MSE reult : " + str(MSE_A))

# 3-(b) Hitung MSE
MSE_B = 0
MSE_B = np.square(np.subtract(dct_blockList, b_idct)).mean(axis=None)
print("3-(b) MSE reult : " + str(MSE_B))

# 3-(c) Hitung MSE
MSE_C = 0
MSE_C = np.square(np.subtract(dct_blockList, c_idct)).mean(axis=None)
print("3-(c) MSE reult : " + str(MSE_C))

# 3-(d) Hitung MSE
MSE_D = 0
MSE_D = np.square(np.subtract(dct_blockList, d_idct)).mean(axis=None)
print("3-(d) MSE reult : " + str(MSE_D))

# 3-(e) Hitung MSE
MSE_E = 0
MSE_E = np.square(np.subtract(dct_blockList, e_idct)).mean(axis=None)
print("3-(e) MSE reult : " + str(MSE_E))

# 3-(a) Perhitungan PSNR
PSNR_A = 0
PSNR_A = 10 * math.log10(255**2 / MSE_A)
print("3-(A) PSNR result : " + str(PSNR_A) + "(dB)")

# 3-(b) Perhitungan PSNR
PSNR_B = 0
PSNR_B = 10 * math.log10(255**2 / MSE_B)
print("3-(B) PSNR result : " + str(PSNR_B) + "(dB)")

# 3-(c) Perhitungan PSNR
PSNR_C = 0
PSNR_C = 10 * math.log10(255**2 / MSE_C)
print("3-(C) PSNR result : " + str(PSNR_C) + "(dB)")

# 3-(d) Perhitungan PSNR
PSNR_D = 0
PSNR_D = 10 * math.log10(255**2 / MSE_D)
print("3-(D) PSNR result : " + str(PSNR_D) + "(dB)")

# 3-(e) Perhitungan PSNR
PSNR_E = 0
PSNR_E = 10 * math.log10(255**2 / MSE_E)
print("3-(E) PSNR result : " + str(PSNR_E) + "(dB)")


# 3-(a) Perhitungan SNR
SNR_A = 0
SNR_A = 10 * math.log10(8**2 / MSE_A)
print("3-(A) SNR result : " + str(SNR_A) + "(dB)")

# 3-(b) Perhitungan SNR
SNR_B = 0
SNR_B = 10 * math.log10(8**2 / MSE_B)
print("3-(B) SNR result : " + str(SNR_B) + "(dB)")

# 3-(c) Perhitungan SNR
SNR_C = 0
SNR_C = 10 * math.log10(8**2 / MSE_C)
print("3-(C) SNR result : " + str(SNR_C) + "(dB)")

# 3-(d) Perhitungan SNR
SNR_D = 0
SNR_D = 10 * math.log10(8**2 / MSE_D)
print("3-(D) SNR result : " + str(SNR_D) + "(dB)")

# 3-(e) Perhitungan SNR
SNR_E = 0
SNR_E = 10 * math.log10(8**2 / MSE_E)
print("3-(E) SNR result : " + str(SNR_E) + "(dB)")


def runLength2bytes(code):
    # perhitungan data input untuk mencari rasio gambar kompresi
    return bytes([len(code) % 8]+[int(code[i:i+8], 2) for i in range(0, len(code), 8)])


# 3-a Hasil kompresi sebelum dan sesudah
code = ""
gambar = imageio.imread('original.gif')
gambar2 = imageio.imread('3(a)_reconstruct_image.jpg')
print("Rasio Kompresi 3 - a dengan gambar asli:                    %.2f : 1" %
      (gambar.size/2**20/(gambar2.size/2**20)))  # perhitungan kompresi gambar

# 3-b Hasil kompresi sebelum dan sesudah
code = ""
gambar3 = imageio.imread('3(a)_reconstruct_image.jpg')
gambar4 = imageio.imread('3(b)_reconstruct_image.jpg')
print("Rasio Kompresi 3 - b dengan 3 - a :                         %.2f : 1" %
      (gambar3.size/2**20/(gambar4.size/2**20)))  # perhitungan kompresi gambar

# 3-c Hasil kompresi sebelum dan sesudah
code = ""
gambar5 = imageio.imread('3(b)_reconstruct_image.jpg')
gambar6 = imageio.imread('3(c)_reconstruct_image.jpg')
print("Rasio Kompresi 3 - c dengan 3 - b :                         %.2f : 1" %
      (gambar5.size/2**20/(gambar6.size/2**20)))  # perhitungan kompresi gambar

# 3-d Hasil kompresi sebelum dan sesudah
code = ""
gambar7 = imageio.imread('3(c)_reconstruct_image.jpg')
gambar8 = imageio.imread('3(d)_reconstruct_image.jpg')
print("Rasio Kompresi 3 - d dengan 3 - c :                         %.2f : 1" %
      (gambar7.size/2**20/(gambar8.size/2**20)))  # perhitungan kompresi gambar

# 3-e Hasil kompresi sebelum dan sesudah
code = ""
gambar9 = imageio.imread('3(d)_reconstruct_image.jpg')
gambar10 = imageio.imread('3(e)_reconstruct_image.jpg')
print("Rasio Kompresi 3 - e dengan 3 - d :                         %.2f : 1" %
      (gambar9.size/2**20/(gambar10.size/2**20)))  # perhitungan kompresi gambar

# 3-e Hasil kompresi sebelum dan sesudah
code = ""
gambar11 = imageio.imread('original.gif')
gambar12 = imageio.imread('3(e)_reconstruct_image.jpg')
print("Rasio Kompresi akhir dengan gambar asli :                   %.2f : 1" %
      (gambar11.size/2**20/(gambar12.size/2**20)))  # perhitungan kompresi gambar