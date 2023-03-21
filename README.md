# BlackScholes
All the BlackScholes Formula Useful to calculate Price/Greeks and Implied/Historical Volatility

Simply call all funcs in the 'bsformula.py' to generate the value you want. follow the input order.

For instance:

```
try_ = bsformula.BlackScholes(500, 500, 0.1, 0.02, 2, 0) 
# inputs  areStock / Strike / Sigma / Risk Free Rate / Time Duration / Dividend
```
```
try_.black_scholes()
print(try_.Call)
print(try_.Put)
try_.delta()
try_.gamma()
try_.vega()
try_.theta()
```

Also you can call other functions similarly as above.

