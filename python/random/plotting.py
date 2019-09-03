

def plot_two_datasets(P1, P2):
    data = []
    col = 0
    for line in range(P1.shape[1]):
        col=col+1
        trace = dict(
            type='scatter',
            x=[P2.transpose()[line, 0], P1.transpose()[line, 0]],
            y=[P2.transpose()[line, 1], P1.transpose()[line, 1]],
            mode='lines+markers',
            marker=dict(
                color=[1,2],
                size=12,
                line=dict(
                    color='rgba(217, 217, 217, 0.14)',
                    width=0.5),
                opacity=0.8)
        )
        data.append(trace)

    layout = dict(
        showlegend=True,
        scene=dict(
            xaxis=dict(title='PC1'),
            yaxis=dict(title='PC2')
        )
    )

    fig = dict(data=data, layout=layout)
    py.offline.plot(fig, filename='test')



def plot_one_dataset(P1, plotname):
    colors = {'Iris-setosa': '#0D76BF',
              'Iris-versicolor': '#00cc96',
              'Iris-virginica': '#EF553B'}
    data=[]
    for name, col in zip(('Iris-setosa', 'Iris-versicolor', 'Iris-virginica'), colors.values()):
        trace = dict(
            type='scatter',
            x=P1.transpose()[labels == name, 0],
            y=P1.transpose()[labels == name, 1],
            mode='markers',
            name=name,
            marker=dict(
                color=col,
                size=12,
                line=dict(
                    color='rgba(217, 217, 217, 0.14)',
                    width=0.5),
                opacity=0.8)
        )
        data.append(trace)

    layout = dict(
        showlegend=True,
        scene=dict(
            xaxis=dict(title='PC1'),
            yaxis=dict(title='PC2')
        )
    )

    fig = dict(data=data, layout=layout)
    py.offline.plot(fig, filename=plotname)
