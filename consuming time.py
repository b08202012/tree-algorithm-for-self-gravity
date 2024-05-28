method = 'no'

if method == 'no':
    import time

    def aa():
        i = 0
        while i<5:
            i = i + 1
            time.sleep(0.5)
            print('A:', i)

    def bb():
        i = 0
        while i<50:
            i = i + 10
            time.sleep(0.5)
            print('B:', i)

    initial_t =  time.time()
    aa()    # 先執行 aa()
    bb()    # aa() 結束後才會執行 bb()
    final_t =  time.time()
    # print("final t", initial_t)
    print('end')
    total_t = final_t - initial_t
    print("total t", total_t, '(s)')
else:

    import threading
    import time

    def aa():
        i = 0
        while i<5:
            i = i + 1
            time.sleep(0.5)
            print('A:', i)

    def bb():
        i = 0
        while i<50:
            i = i + 10
            time.sleep(0.5)
            print('B:', i)

    a = threading.Thread(target=aa)  # 建立新的執行緒
    b = threading.Thread(target=bb)  # 建立新的執行緒

    initial_t =  time.time()
    a.start()  # 啟用執行緒
    b.start()  # 啟用執行緒
    a.join()   # 等到執行完畢，繼續
    b.join()   # 等到執行完畢，繼續
    final_t =  time.time()
    # print("final t", initial_t)
    print('end')
    total_t = final_t - initial_t
    print("total t", total_t, '(s)')