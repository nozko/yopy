#/usr/bin/env python
# coding: utf-8

'''
snow and road-heating simulation

引用論文
(1) main ロードヒーティングの期間融雪負荷シミュレーション
(2)基礎杭利用による地熱融雪法の設計施工運転と数値シミュレーション
(3)広域に適用可能な融雪･積雪水量モデル
'''

from __future__ import print_function
import math
import commands
from datetime import datetime
import time
from string import *
import argparse
import random
import linecache
import sys

### constant ###
pi = 3.14		# 				円周率
temp_s = 0		# [℃ ]			雪氷の融解温度
R = 52			# [W/m^2]		夜間放射量
alpha_r = 5.2	# [W/(m^2*K)]	放射熱伝導率
alpha_m = 233	# [W/(m^2*K)]	雪-路面熱伝導率
a = 0.2			# 				路面または雪面の日射吸収率
epsilon = 0.9	# 				路面または雪面の長波長放射率
#C = 1950		# [kJ/(m^3*K)]	アスファルト比熱
kr = 0.7		# [W/(m*K)]		路盤(アスファルト)の熱伝導率
kp = 0.38		# [W/(m*K)]		パイプの熱伝導率
ks = 0.093		# [W/(m*K)]		雪の熱伝導率
ri = 0.0175		# [m]			パイプ内径
ro = 0.022		# [m]			パイプ外径
cw = 4184		# [J/(kg*K)]	水の比熱
Le = 2267000	# [J/kg]		蒸発潜熱
Lm = 333500		# [J/kg]		融解潜熱
ro_i = 917		# [kg/m^3]		氷の密度
ro_w = 1000		# [kg/m^3]		水の密度
vw = 10			# [m/sec]		温水流速
Qp = 300		# [W/m^2]		供給熱量(パイプ)


### variable ###
#A			パイプの断面積   使ってない

### input ###
#F			降雪または降雨強度
#temp_o		外気温
#RH			外気の湿度
#I			日射量
#V			風速
#P			気圧[Pa]

#on_t		onになってからの経過時間(min)
#off_t		offになってからの経過時間(min)


class simulate:

	def __init__(self):
		weatherF = open('weather.csv', 'r')
		self.weathers = weatherF.readlines()


	# 外気の絶対湿度(Xo)
	def calc_Xo(self, RH, temp_o, P):
		e = 6.11 * ( 10**( 7.5*temp_o/(237.3+temp_o) ) )
		ep = e * RH / 100
		Xo = 0.622*ep / ( P - 0.378*ep )
		print('Xo :', Xo)
		return Xo


	# 路面温度に対する飽和絶対湿度(Xr)
	def calc_Xr(self, temp_o, P):
		es = math.exp(6.414672 + 0.07266115*temp_o)
		Xr = (0.622 * es) / P
		print('Xr :', Xr)
		return Xr


	# 雪表温度に対する飽和絶対湿度(Xs)
	def calc_Xs(self, temp_o, P):
		es = math.exp(6.414548 + 0.08238957*temp_o)
		Xs = (0.622 * es) / P
		print('Xs :', Xs)
		return Xs


	# 降雪の密度(ro_f)
	def calc_ro_f(self, temp_o):
		ro_f = 1000 / (0.091*temp_o*temp_o - 1.81*temp_o + 9.47)
		print('ro_f :', ro_f)
		return ro_f


	# 対流熱伝達率(alpha_c)
	def calc_alpha_c(self, V):
		alpha_c = 7.4 + V
		print('alpha_c :', alpha_c)
		return alpha_c

	# alpha_x
	def calc_alpha_x(self, V):
		alpha_x = 0.0152*V + 0.0178
		print('alpha_x :', alpha_x)
		return alpha_x

	# alpha_i
	def calc_alpha_i(self):
		alpha_i = 5000 + (2/3)*vw
		print('alpha_i :', alpha_i)
		return alpha_i
	

	# 降雪の固相(氷)の割合(rs)
	def calc_rs(self, temp_o):
		if(temp_o > 2):	rs = 0
		elif(temp_o<0):	rs = 1
		else:			rs = 1 - (temp_o/2)
		print('rs :', rs)
		return rs



	# 温水温度(パイプ中間)
	def calc_temp_w(self, on_t):
		temp_w = 0.05 * on_t
		print('temp_w :', temp_w)
		return temp_w


	# 相当外気温(temp_e)
	def calc_temp_e(self, temp_o, I, alpha_c):
		temp_e = temp_o + (a*I - epsilon*R) / (alpha_c + alpha_r)
		print('temp_e :', temp_e)
		return temp_e


	# 融雪量(M)
	def calc_M(self, S_1, elapsed_t, rs, F):
		if(elapsed_t==0):	M = 0
		else:				M = (S_1/elapsed_t) + (rs*F)
		return M

	def recalc_M(self, S_1, elapsed_t, rs, F):
		bigM = False
		if(S_1!=0 and S==0):	bigM = True
		if(W_1!=0 and W==0):	bigM = True
		if(bigM==False and temp_r<=0 and temp_e<=0 and W_1==0):
			M = 0
		return M



	# 積雪の乾燥密度(ro_s)
	def calc_ro_s(self, temp_o, on_t, S_1, F, M, elapsed_t, ro_s_1, ro_f):
		if(temp_o>0 or on_t>0):	#融雪時
			ro_s = (S_1 + F*elapsed_t) / ( (S_1/ro_s_1) + F*elapsed_t/ro_f )
		else:					#凍結時
			if( (S_1/ro_s_1 + F*elapsed_t/ro_f - M*elapsed_t/ro_i)==0 ):
				ro_s = 1
			else:
				ro_s = (S_1 + F*elapsed_t - M*elapsed_t) / ( S_1/ro_s_1 + F*elapsed_t/ro_f - M*elapsed_t/ro_i )
		print('ro_s :', ro_s)
		return ro_s


	# 路面の雪氷量(S)
	def calc_S(self, S_1, rs, F, M, elapsed_t):
		S = S_1 + (rs*F - M)*elapsed_t
		return S


	# 水分の路面からの浸透高さ(dw)
	def calc_dw(self):
		#dw = W / (ro_w - ro_s)
		dw = 0
		print('dw :', dw)
		return dw


	# 路面の積雪深(ds)
	def calc_ds(self, S, ro_s):
		if(ro_s==0):	ds = 0
		else:			ds = S / ro_s
		return ds


	# 路面からの蒸発量(Er)
	def calc_Er(self, alpha_x, Xr, Xo, S):
		if(S>0):	Er = 0
		else:		Er = alpha_x * (Xr - Xo)
		print('Er :', Er)
		return Er


	# 路面から外界への潜熱伝熱量(Qer)
	def calc_Qer(self, Er):
		Qer = Le * Er
		print('Qer :', Qer)
		return Qer


	# 雪から外界への顕熱伝熱量(Qas)
	def calc_Qas(self, alpha_c, ds, dw, temp_e):
		Ks = 1 / ( (1/alpha_c) + (ds-dw)/ks )
		Qas = -1 * Ks * temp_e
		print('Qas :', Qas)
		return Qas
		#alpha_c -> alpha_o ???


	# 雪面からの蒸発量(Es)
	def calc_Es(self, alpha_x, Xs, Xo):
		Es = alpha_x * (Xs - Xo)
		print('Es :', Es)
		return Es

	# 雪から外界への潜熱伝熱量(Qes)
	def calc_Qes(self, Es):
		Qes = Le * Es
		print('Qes :', Qes)
		return Qes

	# 融雪熱量(Qm)
	def calc_Qm(self, M):
		Qm = Lm * M
		return Qm


	# 路面から相変化中の雪への伝熱量(Qs)
	def calc_Qs(self, Qas, Qes, Qm):
		#Qs = alpha_m * (temp_r - temp_s)
		Qs = Qas + Qes + Qm
		print('Qs :', Qs)
		return Qs


	# 路面の温度(temp_r)
	def calc_temp_r(self, Qs, alpha_m):
		'''
		if(on_t>0):
			temp_r = temp_r_1 + (1/30)*on_t
		elif(off_t>0):
			temp_r = temp_r_1 - (1/45)*off_t
		else:
			temp_r = temp_r_1
		'''
		temp_r = Qs/alpha_m + temp_s
		return temp_r


	# 路面から外界への顕熱伝熱量(Qar)
	def calc_Qar(self, temp_r, temp_e, alpha_c, ds):
		Ks = 1 / ( 1/alpha_c + ds/ks )
		Qar = Ks * ( temp_r - temp_e )
		print('Qar :', Qar)
		return Qar
		#alpha_c -> alpha_o???


	# 路面のうち雪氷あるいは水分の相変化に寄与する面積の割合(f)
	def calc_f(self, M, W_1, elapsed_t, rs, F, temp_r, temp_e, Er, Es):
		if(temp_r<=0 and temp_e<=0 and W_1==0):
			f = 0.0
			E = Er
		elif(temp_e>0):
			f = 1.0
			E = Es
		elif(temp_e<=0 and temp_r>0):
			f = 1.0
			E = Es
		elif(temp_e<=0 and temp_r<=0 and W_1>0):
			f = 1.0
			E = Es
		elif(bigM==True):
			E = (W_1/elapsed_t) + (1-rs)*F + M
			f = (E-Er) / (Es-Er)
		print('f :', f)
		print('E :', E)
		return f, E


	# 路面放熱量(Qr)
	def calc_Qr(self, f, Qs, Qar, Qer):
		Qr = f*Qs + ( 1.0 - f )*( Qar + Qer )
		print('Qr :', Qr)
		return Qr


	# 路面の水分量(W)
	def calc_W(self, W_1, rs, F, M, E, elapsed_t):
		W = W_1 + ( (1-rs)*F + M - E )*elapsed_t
		return W


	# パイプ放熱量(Qp)
	def calc_Qp(self, temp_w, temp_r, alpha_i):
		Kp = 1 / ( 1/(2*pi*alpha_i*ri) + (math.log(ro/ri)/(2*pi*kp)) )
		Qp = Kp * ( temp_w - temp_r )
		print('Qp :', Qp)
		return Qp



	# 天気情報
	def get_weather(self, num):
		w = self.weathers[num]
		weather = w.split(', ')
		return weather


	# 結果出力
	def result_ouput(self, datetime, M, S, ds, Qm, temp_r, W):
		outf = open('result_simulate.txt', 'a')
		result = str(datetime)+', '+str(M)+', '+str(S)+', '+str(ds)+', '+str(Qm)+', '+str(temp_r)+', '+str(W)+'\n'
		print(result)
		outf.write(result)



if __name__ == '__main__':

	num = 0
	weather = simulate().get_weather(num)
	print(weather)

	## initial states ##
	datetime = datetime.strptime(weather[0], "%m-%d %H:%M")
	F = float(weather[1])					#降雪または降雨強度
	temp_o = float(weather[2])				#外気温
	RH = float(weather[3])					#外気の湿度
	I = float(weather[4])					#日射量
	V = float(weather[5])					#風速
	P = float(weather[6].split('\n')[0])	#気圧[Pa]

	on_t = 0
	off_t = 0

	Xo = simulate().calc_Xo(RH, temp_o, P)


	Xr = simulate().calc_Xr(temp_o, P)
	Xs = simulate().calc_Xs(temp_o, P)
	ro_f = simulate().calc_ro_f(temp_o)
	alpha_c = simulate().calc_alpha_c(V)
	alpha_x = simulate().calc_alpha_x(V)
	alpha_i = simulate().calc_alpha_i()
	rs = simulate().calc_rs(temp_o)
	temp_w = simulate().calc_temp_w(on_t)

	temp_e = simulate().calc_temp_e(temp_o, I, alpha_c)
	M = 0
	ro_s = simulate().calc_ro_s(temp_o, on_t, 0, F, M, 0, 1, ro_f)
	S = 0
	dw = simulate().calc_dw()
	ds = simulate().calc_ds(S, ro_s)
	Er = simulate().calc_Er(alpha_x, Xr, Xo, S)
	Qer = simulate().calc_Qer(Er)
	Qas = simulate().calc_Qas(alpha_c, ds, dw, temp_e)
	Es = simulate().calc_Es(alpha_x, Xs, Xo)
	Qes = simulate().calc_Qes(Es)
	Qm = simulate().calc_Qm(M)
	Qs = simulate().calc_Qs(Qas, Qes, Qm)
	temp_r = simulate().calc_temp_r(alpha_m, Qs)
	Qar = simulate().calc_Qar(temp_r, temp_e, alpha_c, ds)
	f, E = simulate().calc_f(M, 0, 0, rs, F, temp_r, temp_e, Er, Es)
	Qr = simulate().calc_Qr(f, Qs, Qar, Qer)
	W = 0
#	Qp = simulate().calc_Qp(temp_w, temp_r, alpha_i)

	b_time = datetime
	S_1 = S
	W_1 = W
	ro_s_1 = ro_s

	simulate().result_ouput(datetime, M, S, ds, Qm, temp_r, W)

	roopnum = len(simulate().weathers)
	for n in range(roopnum-1):

		time.sleep(5)

		num += 1
		weather = simulate().get_weather(num)
		print(weather)
	
		## initial states ##
		datetime = datetime.strptime(weather[0], "%m-%d %H:%M")
		F = float(weather[1])					#降雪または降雨強度
		temp_o = float(weather[2])				#外気温
		RH = float(weather[3])					#外気の湿度
		I = float(weather[4])					#日射量
		V = float(weather[5])					#風速
		P = float(weather[6].split('\n')[0])	#気圧[Pa]
	
		on_t = 0
		off_t = 0
	
		Xo = simulate().calc_Xo(RH, temp_o, P)
		elapsed_time = datetime - b_time
		print('elapsed time :', elapsed_time)
		elapsed_t = (elapsed_time.seconds)/60
		print('elapsed minutes :', elapsed_t)
		

		Xr = simulate().calc_Xr(temp_o, P)
		Xs = simulate().calc_Xs(temp_o, P)
		ro_f = simulate().calc_ro_f(temp_o)
		alpha_c = simulate().calc_alpha_c(V)
		alpha_x = simulate().calc_alpha_x(V)
		alpha_i = simulate().calc_alpha_i()
		rs = simulate().calc_rs(temp_o)
		temp_w = simulate().calc_temp_w(on_t)
	
		temp_e = simulate().calc_temp_e(temp_o, I, alpha_c)
		M = simulate().calc_M(S_1, elapsed_t, rs, F)
		ro_s = simulate().calc_ro_s(temp_o, on_t, S_1, F, M, elapsed_t, ro_s_1, ro_f)
		S = simulate().calc_S(S_1, rs, F, M, elapsed_t)
		dw = simulate().calc_dw()
		ds = simulate().calc_ds(S, ro_s)
		Er = simulate().calc_Er(alpha_x, Xr, Xo, S)
		Qer = simulate().calc_Qer(Er)
		Qas = simulate().calc_Qas(alpha_c, ds, dw, temp_e)
		Es = simulate().calc_Es(alpha_x, Xs, Xo)
		Qes = simulate().calc_Qes(Es)
		Qm = simulate().calc_Qm(M)
		Qs = simulate().calc_Qs(Qas, Qes, Qm)
		temp_r = simulate().calc_temp_r(alpha_m, Qs)
		Qar = simulate().calc_Qar(temp_r, temp_e, alpha_c, ds)
		f, E = simulate().calc_f(M, W_1, elapsed_t, rs, F, temp_r, temp_e, Er, Es)
		Qr = simulate().calc_Qr(f, Qs, Qar, Qer)
		W = simulate().calc_W(W_1, rs, F, M, E, elapsed_t)
#		Qp = simulate().calc_Qp(temp_w, temp_r, alpha_i)

		b_time = datetime
		S_1 = S
		W_1 = W
		ro_s_1 = ro_s
		
		simulate().result_ouput(datetime, M, S, ds, Qm, temp_r, W)
