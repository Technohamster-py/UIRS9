import geocoder


city = "Novosibirsk"

if __name__ == '__main__':
    city_cords = geocoder.google(city).latlng