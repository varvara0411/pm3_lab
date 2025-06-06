import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('times.csv')

plt.figure()
plt.plot(df['N'], df['LCG'], label='LCG')
plt.plot(df['N'], df['Xorshift'], label='Xorshift32')
plt.plot(df['N'], df['MiddleSquare'], label='MiddleSquare')
plt.plot(df['N'], df['MT19937'], label='std::mt19937')
plt.xlabel('Размер выборки (N)')
plt.ylabel('Время генерации (микросекунды)')
plt.title('Зависимость времени генерации от размера выборки')
plt.legend()
plt.tight_layout()
plt.show()