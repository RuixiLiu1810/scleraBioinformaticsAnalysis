import requests

with open("1.mp4",'wb') as f:
    hd = {
        "User-Agent": "Mozilla/5.0 (X11; Linux aarch64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/144.0.0.0 Safari/537.36 CrKey/1.54.250320",
        "Referer": "https://meeting.tencent.com/"
    }
    url = "https://ylz.cos.meeting.tencent.com/cos/200000001/2014588190486683648/2014588190486683649/TM-20260123141950-406142634-recording-1.mp4?token=eJxskElz2jAYhv-Lr9RYyydZZqYHnDBkMQmLocCFUbwq4N0qMaX_vUMCnR6q46tve55fhu8t-rIsd1qFxsCQeZF3mWo7E4JmfL8N110RV27kRnW3fSEJr_RHuHls-dL49tUatKrIjYExPcjunLZt2Qwsq9P5QZ9SZcqTiQkTFIA7Tj9Wh6ifdVVwKHTYD4rMCorGIujrYYsgDEwI7CAQnAvKQfwncyx_YhJEOMKEYsAOQyYgjoFwCmYdBUUdqjwxcT8r4YzO-IyckHF4Iya8gW0CIGoKzphJA8EEAhlJyq88TbjfybL8lIEFALcBkE1t5DByWXAta7syMgbGfHT3Or9_fBnfYpVdYmxzh1HMibhNVYkxMO6CtXtKm656ou66Q_NtufCeD5MA3-tC8sSOH9xqH3vNdto2uO6xk6V1b_byvJ8qfxOx-KdXHkt__DB8SMnEqWk8o0td93r28TEnoTWeeNPFMJTPI49vF-kx3qz5GFXJaqXnvu1XCyD5bOVwXcQyO70u7R8Wx40LSzaU7-J9A-EoV6NCk6f58Pv18OijVHW0k3Eb1TcwRslfsFS1u6SW3eXvH4Of-ozffwIAAP__EG6yuw"
    rq = requests.get(url=url, headers=hd)
    f.write(rq.content)