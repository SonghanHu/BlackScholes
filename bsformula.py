import numpy as np
from scipy import stats


class BlackScholes:

    def __init__(self, s, e, sigma, r, delta_t, dividend=0):

        self.S = s
        self.E = e
        self.sigma = sigma
        self.r = r
        self.delta_t = delta_t
        self.dividend = dividend
        partial = self.r - self.dividend + (1/2)*self.sigma**2
        self.d1 = (np.log(self.S/self.E)+partial*self.delta_t)/(self.sigma*np.sqrt(self.delta_t))
        self.d2 = self.d1 - self.sigma*np.sqrt(self.delta_t)

    def black_scholes(self):

        Part_E = self.E * np.exp(-self.r * self.delta_t) * stats.norm.cdf(self.d2)
        Part_S = self.S * np.exp(-self.dividend * self.delta_t) * stats.norm.cdf(self.d1)
        self.Call = Part_S - Part_E
        self.Put = self.Call - self.S + self.E * np.exp(-self.r * self.delta_t)

    def delta(self):

        self.Cdelta = stats.norm.cdf(self.d1)
        self.Pdelta = -stats.norm.cdf(-self.d1)
        print(f'Call Delta: {self.Cdelta}; Put Delta: {self.Pdelta}')

    def vega(self):

        self.Cvega = np.sqrt(self.delta_t) * self.S * stats.norm.pdf(self.d1)
        self.Pvega = self.Cvega
        print(f'Call Vega: {self.Cvega}; Put Vega: {self.Pvega}')

    def theta(self):

        p = ((self.sigma/np.sqrt(self.delta_t)/2)*stats.norm.pdf(self.d1))*self.S
        self.Ctheta = -self.r*self.E*np.exp(-self.r*self.delta_t)*stats.norm.cdf(self.d2) - p
        self.Ptheta = self.r*self.E*np.exp(-self.r*self.delta_t)*stats.norm.cdf(-self.d2) - p
        print(f'Call Theta: {self.Ctheta}; Put Vega: {self.Ptheta}')

    def gamma(self):

        self.Cgamma = stats.norm.pdf(self.d1) / self.S * self.sigma * np.sqrt(self.delta_t)
        self.Pgamma = self.Cgamma
        print(f'Call Gamma: {self.Cgamma}; Put Gamma: {self.Pgamma}')
        

class implied_vol:

    def __init__(self, S, E, r, delta_t, Call, Put) -> None:

        self.C = Call
        self.P = Put
        self.S = S
        self.E = E
        self.r = r
        self.delta_t = delta_t

    def implied_vol_call(self):
        
        volatility_candidates = np.arange(0.01,4,0.0001)
        price_differences = np.zeros_like(volatility_candidates)

        for i in range(len(volatility_candidates)):
            candidate = volatility_candidates[i]
            a = BlackScholes(self.S,self.E,candidate,self.r,self.delta_t)
            a.black_scholes()
            price_differences[i] = self.C - a.Call
        idx = np.argmin(abs(price_differences))
        self.implied_volatility = volatility_candidates[idx]
        return self.implied_volatility
    
    def implied_vol_put(self):

        volatility_candidates = np.arange(0.01,4,0.0001)
        price_differences = np.zeros_like(volatility_candidates)

        for i in range(len(volatility_candidates)):
            candidate = volatility_candidates[i]
            a = BlackScholes(self.S, self.E, candidate, self.r, self.delta_t)
            a.black_scholes()
            price_differences[i] = self.P - a.Put
        idx = np.argmin(abs(price_differences))
        self.implied_volatility = volatility_candidates[idx]
        return self.implied_volatility

class date_to_day:

    def __init__(self,year,month,day) -> None:

        self.a = year
        self.b = month
        self.c = day

    def function(self):
        if int(self.b) == 1:
            return self.c
        else:
            if (int(self.a)%400 != 0 or int(self.a)%4 != 0):
                month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
                days = int(self.c)
                for i in range(int(self.b)-1):
                    days += month[i]
            else:
                month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
                days = int(self.c)
                for i in range(int(self.b)-1):
                    days += month[i]
            return days

class logmode:

    def __init__(self, s, mu, sigma, delta_t) -> None:

        self.S = s
        self.mu = mu
        self.sigma = sigma
        self.delta_t = delta_t

    def logmode(self):

        self.mode = self.S * np.exp((self.mu - (3/2)*self.sigma**2)*self.delta_t)
        self.max = (1/(np.sqrt(2*np.pi)*self.S*self.sigma*np.sqrt(self.delta_t)))*np.exp((-self.mu+self.sigma**2)*(self.delta_t))
        print(f'The stock price is most likely to be {self.mode} with probability{self.max}')

    def proba(self, s):

        self.mu1 = np.log(self.S) + (self.mu - (1/2) *self.sigma**2)*self.delta_t
        self.sigma1 = self.sigma * np.sqrt(self.delta_t)
        x = (np.log(s) - self.mu1)/self.sigma1
        self.proba_less = stats.norm.cdf(x)
        self.proba_more = 1 - self.proba_less
        # print(f"Call_Out/Put_In:P = {self.proba_less},Call_In/ Put_Out: P = {self.proba_more}")

class historical_val:

    def __init__(self, price_list) -> None:

        self.price_list = price_list
        self.return_list = []

    def calculation(self):

        for i in range(1, len(self.price_list)):
            a = np.log(((self.price_list[i]-self.price_list[i-1])/self.price_list[i-1] + 1))
            self.return_list.append(a)
        self.return_array = self.return_list
        self.mu = np.mean(self.return_array)
        self.sigma = np.sqrt((1 / (len(self.return_list) - 1)) * np.sum(((self.return_array - self.mu) ** 2)))
        self.mu_year = self.mu * 252
        self.sigma_year = self.sigma * np.sqrt(252)

class bet:

    def __init__(self):

        pass

    def calc_p(self,Oa,Ob):
        self.x_1 = (1 + Ob ) / (2 + Oa + Ob)
        self.p = (Oa * Ob-1)/(2 + Oa + Ob)
        self.win_a = Oa/(1+Oa)
        self.win_b = Ob/(1+Ob)
        return (self.x_1, self.p)

    def calc_O(self, win_a,win_b):

        self.Oa = win_a/(1-win_a)
        self.Ob = win_b/(1-win_b)

