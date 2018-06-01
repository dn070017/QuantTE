 # QuantTE
 ## About
 QuantTE is an analsis pipeline that utilize the <b>RNA-Seq</b> data and <b>RepeatMasker table</b> to quantify the repeat element. There are two main stages for QuantTE: <b>extraction</b> and <b>quantification</b>. The first stage used the RepeatMasker table and genome sequences to extract the transposable element in interest while the second stage used Kallisto to quantify the abundance of transposable element.
 ## Manual
 - Extract stage:
 ```shell
 python3 QuantTE/main.py --extract --input input.json
 ```
 - Quant stage:
 ```shell
 python3 QuantTE/main.py --quant --input input.json
 ```
 - Run Complete Analysis
  ```shell
 python3 QuantTE/main.py --extract --quant --input input.json
 ```
