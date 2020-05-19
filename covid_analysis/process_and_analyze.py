import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime

class process_and_analyze:
    
    def __init__(self, time_range = -1):
        self.dates = self.get_dates(time_range)
        
    def get_dates(self, time_range = -1):
        if time_range == -1:
            return pd.date_range(start = '1/22/2020', end = pd.datetime.today()-datetime.timedelta(days = 1)).strftime("%m-%d-%Y").to_list()
        else:
            return pd.date_range(start = pd.datetime.today()-datetime.timedelta(days = time_range + 1), end = pd.datetime.today()-datetime.timedelta(days = 1)).strftime("%m-%d-%Y").to_list()

    def get_raw_data(self):
        dates = self.dates
        dataframes = {}
        for date in dates:
            try:
                dataframes[date] = pd.read_csv('/Users/vahsengrouplaptop/covid/COVID-19/csse_covid_19_data/csse_covid_19_daily_reports/%s.csv'%(date))
            except FileNotFoundError:
                continue
        
            dataframes[date] = dataframes[date].rename(columns={"Country_Region": "Country/Region", "Province_State": "Province/State"}) #keep format the same for all dates
        return dataframes

    def get_names_of_states(self):
        states = {'Alabama': 'AL','Alaska': 'AK','Arizona': 'AZ','Arkansas': 'AR','California': 'CA','Colorado': 'CO','Connecticut': 'CT','Delaware': 'DE','District of Columbia': 'DC','Florida': 'FL','Georgia': 'GA','Guam': 'GU','Hawaii': 'HI','Idaho': 'ID','Illinois': 'IL','Indiana': 'IN','Iowa': 'IA','Kansas': 'KS', 'Kentucky': 'KY', 'Louisiana': 'LA','Maine': 'ME','Maryland': 'MD','Massachusetts': 'MA','Michigan': 'MI','Minnesota': 'MN','Mississippi': 'MS','Missouri': 'MO','Montana': 'MT','Nebraska': 'NE','Nevada': 'NV','New Hampshire': 'NH','New Jersey': 'NJ','New Mexico': 'NM','New York': 'NY','North Carolina': 'NC','North Dakota': 'ND','Northern Mariana Islands':'MP','Ohio': 'OH','Oklahoma': 'OK','Oregon': 'OR','Pennsylvania': 'PA','Puerto Rico': 'PR','Rhode Island': 'RI','South Carolina': 'SC','South Dakota': 'SD','Tennessee': 'TN','Texas': 'TX','Utah': 'UT','Vermont': 'VT','Virginia': 'VA','Washington': 'WA','West Virginia': 'WV','Wyoming': 'WY'}
        states = dict(map(reversed, states.items()))
        df = pd.DataFrame()
        df['abbr'] = list(states.keys())
        df['name'] = [states[key] for key in states.keys()]
        return df

    def get_observables(self, country, *args):
        if len(args)==1:
            data = self.get_raw_data()
            col = 'Province/State'
            val =  args[0]
            states = self.get_names_of_states()
            states = states.loc[(states['name'] == val) | (states['abbr'] == val)]
            abbr, state = [val for val in states['abbr']][0], [val for val in states['name']][0]
            val = state
            for date in data.keys():
                data[date]=data[date].replace(to_replace='.*, %s'%(abbr), value='%s'%(state), regex=True)
        elif len(args) == 0:
            data = self.get_raw_data()
            col = 'Country/Region'
            val = country
            if country == 'China':
                for date in data.keys():
                    data[date]=data[date].replace(to_replace="Mainland China", value="China", regex=False)
        else:
            print("Please enter as arguments either (country) or ('US', state)")
        if country == "World":
            confirmed_case = [data[date]['Confirmed'].sum() for date in data.keys() if data[date]['Confirmed'].sum() > 0]
            death = [data[date]['Deaths'].sum() for date in data.keys() if data[date]['Deaths'].sum() > 0]
        else:
            confirmed_case = [data[date].loc[data[date][col]==val]['Confirmed'].sum() for date in data.keys() if data[date].loc[data[date][col]==val]['Confirmed'].sum() > 0]
            death = [data[date].loc[data[date][col]==val]['Deaths'].sum() for date in data.keys() if data[date].loc[data[date][col]==val]['Deaths'].sum() > 0]
        confirmed_cases = [np.nan for i in range(0,int(len(data.keys()) - len(confirmed_case)))] + confirmed_case
        deaths = [np.nan for i in range(0,int(len(data.keys()) - len(death)))] + death
        observables = pd.DataFrame()
        observables['date'] = data.keys()
        observables['Confirmed'] = confirmed_cases
        observables['Deaths'] = deaths
        observables['Confirmed/day'] = observables['Confirmed'].diff()
        observables['Deaths/day'] = observables['Deaths'].diff()
        return observables

    def fit_data(self, col, country, *args):
        data = self.get_observables(country, *args)
        if (col != 'Confirmed') and (col != 'Deaths') and (col != 'Confirmed/day') and (col != 'Deaths/day'):
            print("Enter 'Confirmed', 'Deaths', 'Confirmed/day', or 'Deaths/day' as the second argument")
        data = data[data[col].notna()]
        x = np.linspace(0,len(data)-1,len(data))
        fit = np.polyfit(x, np.log(data[col].to_numpy()), 1)
        return fit

    def plot(self, col, country, *args):
        data = self.get_observables(country, *args)
        data = data.replace(to_replace ='-2020', value='', regex=True)
        data = data[data[col].notna()]
        plt.subplot(2,1,1)
        if len(args) == 1:
            plt.plot(data['date'], data[col], 'o', label = '%s for %s'%(col, args[0]))
        else:
            plt.plot(data['date'], data[col], 'o', label = '%s for %s'%(col, country))
        if country == 'US':
            x = np.linspace(0,len(data)+2,len(data)+3)
            fit = self.fit_data(col, country, *args)
            double = np.log(2)/fit[0]
            plt.plot(x, np.exp(fit[0]*x+fit[1]), '--', color = 'red', label = 'Exp. fit for US, doubling time = %s days'%(float('%.3g' % double)))
        plt.xticks(np.arange(0, len(data['date']), step=4))
        plt.xlabel('date')
        plt.ylabel(col)
        plt.title('COVID-19 %s'%(col))
        plt.legend()
        plt.subplot(2,1,2)
        if len(args) == 1:
            plt.plot(data['date'], data[col], 'o', label = '%s for %s'%(col, args[0]))
        else:
            plt.plot(data['date'], data[col], 'o', label = '%s for %s'%(col, country))
        #x = np.linspace(0,len(data)+2,len(data)+3)
        #fit = self.fit_data(col, country, *args)
        #plt.plot(x, np.exp(fit[0]*x+fit[1]), '--', color = 'red', label = 'Exp. fit for US, doubling time = %s days'%(float('%.3g' % double)))
        plt.xticks(np.arange(0, len(data['date']), step=4))
        plt.xlabel('date')
        plt.ylabel(col)
        plt.yscale("log")
        plt.legend()
        #plt.show()
        
        
        
