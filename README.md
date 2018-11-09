# DSCC
DSCC源代碼，達世鏈源代碼，達世鏈簡介
DSCC 達世鏈柔支付技術（RouPay）白皮書 DSCC—區塊鏈支付領域的 PayPal！ DSCC即乙太國際支付,基於乙太坊公鏈，在EOS基礎上優化,可跨鏈交易，對稱側鏈+閃電網路,跨 交易所，支援各類主流數位貨幣的幣種兌換。達世鏈 團隊目前擁有彩色區塊鏈技術專利與數個軟體著作版權，達世鏈（DSCC） 團隊的核 心成員均來自于區塊鏈業內，擁有深厚的業內資源及背景。  達世鏈（DSCC） 開發團隊在開發上有著諸多技術創新，由 達世鏈（DSCC） 自主研發的柔支 付技術（RouPay）和 MHT 技術（Matching hedge Technology 匹配對沖技術）等都是由 達世 鏈（DSCC） 團隊自主研發的創新功能。  柔支付技術（RouPay）是由達世鏈（DSCC） 開發團隊自主研發，而基於柔支付技術 （RouPay）為底層打造的柔支付網路（RouPay Network），綜合運用了 2-of-2 多重簽名、鎖 定時間交易、交易構造延後廣播等技術，可以在不需信任的情況，實現區塊鏈資產的零手續費秒 速轉移，在速度、安全性和隱私性方面，足以媲美閃電網路（Lightning Network）
Repository for the Dynamically Skewed Compressed Cache

This does not work as I haven't updated with the latest modifications. I have basically finished implementing locally cache compression on ZSim (lacks one thing on writebacks only), but I've abandoned it to implement the same functionality on Gem5. If you are interested on cache compression in ZSim or Gem5, send me a message and I can try to help you out.

TODO modify timingCache.cpp to make it similar to cache.cpp (as of now it probably doesn't work with compressed caches)

To start:

sudo /sbin/sysctl -w kernel.shmmax=1073741824

export PINPATH=~/pin-2.14-71313-gcc.4.4.7-linux

scons -j16

./build/opt/zsim tests/dscc.cfg
